# Fixing the Membrane Force Loop Without Losing Accuracy — What Went Wrong Before, and a Step-by-Step Fix

*Follow-up to `BM_Complexity_Analysis.md` — grounded in `custom_modules/Archive.cpp`, the prior level-set attempt*

## 1. Was the literature pointing the wrong way? No — but it answers a different question than the one you actually have

Your instinct that something was missing is correct, but it isn't that fast marching is the wrong idea in general. It's that **`Archive.cpp` never actually implements fast marching**, and the grid-based approach it does implement has real, fixable bugs that have nothing to do with grid fineness. Worth separating three distinct design choices that got bundled together in the earlier attempt:

1. *How do you represent "distance to membrane everywhere"* — as a grid (signed distance field), or exactly from the polyline on demand.
2. *If you use a grid, how do you build it fast* — brute force (what you have), or a real redistancing algorithm like Fast Marching (Sethian) / Fast Sweeping (Zhao).
3. *How do you turn "distance to membrane" back into a force on one cell* — read the grid, or compute directly against the polyline.

`Archive.cpp` picked "grid" for (1), never got to (2) at all (`rebuild_signed_distance_field()` is still brute force — the 80KB `fast_marching_method.hpp` sitting in your repo is never called), and used a crude version of (3). So the disappointing accuracy you saw wasn't fast marching failing — fast marching was never actually running. It was the grid-representation-plus-crude-sampling combination in (1)+(3) that hurt you. That's fixable, and fixing it points toward a conclusion you'll probably find reassuring: **for the force computation specifically, you don't need a grid or FMM at all.** More on why below.

## 2. What was actually wrong in `Archive.cpp` — four specific, independent bugs

Reading `cell_interactions_LSM`, `basement_membrane_interactions_LSM`, and `voxel_indices`, here's what's happening, in order of how much damage each one does:

**Bug 1 — the cell's force was evaluated at the wrong point entirely.** Look at `voxel_indices()`:

```cpp
int v = pCell->get_current_voxel_index();
auto& vox = microenvironment.mesh.voxels[v];
double x = vox.center[0];   // <-- the voxel's center, not the cell's position
double y = vox.center[1];
```

This doesn't look up the level-set grid cell containing the cell's actual `(x, y)` position. It looks up which *mechanics voxel* the cell is in, then uses **that voxel's center** as the query point. Your `PhysiCell_settings.xml` has `dx = dy = 20` microns, so a cell can sit up to `√2 × 10 ≈ 14` microns away from the point actually being queried — before you even reach a grid-resolution question. Every distance and every direction fed into the force was computed at a phantom location up to 14 microns from where the cell really is. Against an interaction length of 100 microns that's not catastrophic on its own, but it's compounding with the next three issues, and it's pure unforced error — nothing needs `get_current_voxel_index()` here at all; you already have `pCell->position[0]`, `pCell->position[1]`.

**Bug 2 — nearest-neighbor sampling, not interpolation.** `voxel_indices` does `floor((x - ls_xmin)/ls_dx)` and stops — it snaps to whichever grid cell contains the point and reads `level_set_phi[i][j]` directly, with no bilinear interpolation between neighboring grid values. That makes distance (and therefore force magnitude and direction) **piecewise constant within each 20-micron voxel**, jumping discontinuously every time a cell crosses a voxel boundary. A moving cell would feel a force that's constant, then snaps to a different value, then constant again — exactly the kind of jitter that reads as "the interactions weren't very accurate," independent of how fine the grid is, because the discontinuity is inherent to nearest-neighbor sampling, not to voxel size.

**Bug 3 — the grid spacing was inherited, not chosen.** In `initialize_level_set_duct` (still true in the current `Membrane.cpp`): `ls_dx = mesh.dx; ls_dy = mesh.dy;` — the level-set grid resolution is silently borrowed from BioFVM's diffusion/mechanics mesh, which was sized for solving transport PDEs cheaply, not for resolving membrane geometry. Meanwhile `membrane_num_points = 1000` on a 300-micron-radius circle gives a polyline segment length of about **1.9 microns** — your membrane is already discretized 10× finer than the grid you were about to resample it onto. This is exactly your "I didn't see a solution other than a finer mesh" instinct, and you're right that a uniformly finer mesh would be expensive — a 10× finer grid to match the polyline's own resolution is 100× more voxels in 2-D. That's the real cost blowup you were worried about, and it's avoidable, not inevitable (§3).

**Bug 4 — finite-difference gradients are unreliable exactly where your model gets interesting.** `level_set_gradient` estimates the surface normal via central differences on `φ`. A signed distance function is Lipschitz-continuous but **not differentiable** on the shape's medial axis / skeleton — the set of points equidistant from two or more different boundary features. Concave regions (like the inward dimples your proliferation-driven deformation produces!) push that skeleton closer to the interface, so finite-difference gradients become unreliable precisely in the geometric regime you care most about, and this problem does not go away with a finer grid — it's a structural property of distance functions, not a discretization artifact. This is a well-documented pitfall in the level-set literature, not something specific to your code.

None of these four are arguments against grids or FMM in general. They're arguments against using an *interpolated grid field* as the source of truth for a *force* that needs to be accurate at sub-cell scale. That's the gap in the earlier framework.

## 3. The resolution: stop asking the grid to do the force's job

Here is the reframe that dissolves the dilemma. You have two genuinely different questions that got conflated into one grid:

- **"What force should this specific cell feel right now?"** — needs to be *exact*, because it drives the dynamics you're validating.
- **"Roughly how does distance-to-membrane vary in general?"** — useful for things like phenotype gating, visualization, or a coarse "is this cell even worth checking" pre-filter, where sub-micron accuracy is not required.

Your brute-force `project_point_onto_boundary` already answers the first question *exactly* — that's the code path you trust, and it produces the results you're happy with. The only problem with it is speed, not correctness. So the fix is not "replace exact-but-slow with approximate-but-fast." It's **"keep exact, make it fast by only checking the segments that could possibly matter."** That's a spatial-binning (cell-list) problem, the same one described in the earlier report, and it requires no grid, no interpolation, and no gradient estimation at all — so none of the four bugs above can recur.

A genuine FMM/FSM-built signed distance field is still worth having *for the second question only*, decoupled entirely from the force computation, at whatever resolution and update cadence you like, since nothing dynamically important depends on its precision. That part of the plan doesn't go away — it just stops being on the critical path for "does the model still behave the way I validated it to behave."

## 4. Step-by-step: exact, fast membrane forces via spatial binning

This directly replaces the O(N_cells × N_membrane_nodes) cost in `cell_interactions_cc` with something close to O(N_cells + N_membrane_nodes), while returning bit-identical answers to your current brute-force code (up to floating-point rounding) — meaning validation against your existing trusted results should show effectively zero deviation, not "close enough."

**Step 1 — Pick a bin width equal to your interaction length.** Set the bin width `w = membrane_interaction_length` (currently 100 microns in your config). This is the key correctness guarantee: if `w ≥ L`, then any segment within distance `L` of a query point is guaranteed to lie in the query point's own bin or one of its 8 immediate neighbors (a 3×3 block). No segment farther away can matter, since your forces are already clamped off at `L`. This isn't an approximation — it's the same reasoning PhysiCell's own `agent_grid` relies on for cell-cell mechanics.

**Step 2 — Maintain a bucket list of membrane segments.** Add a simple structure: a 2-D array (or hash map, since the domain is bounded and small) of `vector<int>` — bucket `(i, j)` holds the indices `k` of every membrane segment (defined by nodes `k` and `k+1`) whose position falls in that bucket. Since your segments (≈1.9 microns) are tiny compared to the bin width (100 microns), binning by either the segment's midpoint or its first endpoint is sufficient — you don't need exact segment-vs-box overlap logic.

**Step 3 — Rebuild the buckets once per step, right after the node positions update.** This costs O(N_membrane_nodes) (clear the buckets, loop over all `Np` segments once) — cheap, since it's linear in the membrane's own size, not the product of two populations. Do this at the top of `update_basement_membrane_deformation`, right after `boundary_membrane_pts[i] += node_forces[i] * dt`, before the next step's cell queries need it.

**Step 4 — Replace the linear scan inside `project_point_onto_boundary` with a bucket-restricted scan.** Keep every line of the existing distance/projection math exactly as it is — you're not changing the geometry calculation, only which segments it loops over. Compute the query point's own bucket `(bi, bj)`, gather the segment indices from buckets `(bi±1, bj±1)` (9 buckets total), and run the existing `for (k : candidates)` loop over just that (typically small) candidate list instead of all `Np` segments.

**Step 5 — Replace `is_inside`'s ray cast with a byproduct of the same nearest-segment search.** Right now `is_inside` does its own independent O(Np) ray-casting pass. You don't need it: once you know the nearest segment and projection point from Step 4, the sign of `(query − projection) · outward_normal(best_segment)` tells you inside vs. outside directly — this is in fact the standard way signed distance fields are built in the first place. One-time setup: fix the outward-normal orientation convention by checking it against a known-interior point (e.g., your shape's centroid) once at initialization, so you know whether to rotate each segment's tangent +90° or −90° to get "outward." After that it's a single dot product, replacing another O(Np) loop with an O(1) one. (Keep the old `is_inside` around as an occasional offline sanity check if it'd help you sleep at night — just take it out of the per-step hot path.)

**Step 6 — Swap call sites, change nothing else.** `cell_interactions_cc`, `parallel_cell_division`, and anywhere else that calls `project_point_onto_boundary` or `is_inside` keep the same function signatures and return values — only the internal implementation changes from "scan everything" to "scan the 3×3 neighborhood." No downstream force-law code needs to change at all.

**Step 7 — Validate by exact-match regression, not approximate agreement.** Because this method is exact (not grid-interpolated), the right test is a strict one: run a fixed small scenario (same seed, same few cells, same few membrane nodes) through the current brute-force code and the new binned code, and diff `boundary_membrane_pts` trajectory step by step. Expect agreement to floating-point tolerance (~1e-9), not just "visually similar." If you see drift larger than that, it means a segment near a bucket boundary is being missed — usually a sign that the bin width needs to be a touch larger than `L` (add ~10–20% margin) rather than exactly equal to it, to absorb floating-point edge cases at bucket boundaries.

## 5. Where FMM still fits — later, optional, and no longer risky

Once Steps 1–7 are in and validated, the force loop is no longer the bottleneck driven by grid choices at all, and you're free to revisit a real FMM/FSM-built signed distance field purely for the *other* question (§3) — phenotype-rule distance signals, visualization, or a coarse broad-phase filter — at whatever resolution you like, independent of the mechanics mesh's 20-micron spacing, seeded properly from exact near-interface values, and without any pressure to make it sub-cell-accurate, because nothing dynamically important reads it anymore. That's a much lower-stakes project than trying to make one grid serve both roles at once, which is what made the earlier attempt feel like it was fighting itself.
