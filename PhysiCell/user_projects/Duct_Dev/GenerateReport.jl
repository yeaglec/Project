#=
GenerateReport.jl — PCMM-Integrated Simulation Comparison Report Generator
=============================================================================
Generates an interactive HTML dashboard comparing PhysiCell simulation snapshots,
using PCMM's built-in APIs for parameter detection, cell counting, and data access.

SECTIONS
  Each report is divided into sections. The reliable way to get one section per
  experiment is to pass them explicitly via `sections`:

      GenerateReport("Reports/Report_1_Strain_Force";
          sections = ["Linear"      => out_strain_lin,    # a run()/createTrial() output
                      "Exponential" => out_strain_exp],   # ...or a Vector{Int} of sim IDs
          snapshot_selection = :by_cell_count,
          target_cell_count  = 200)

  For a single, unlabeled section just pass `sim_ids` (or nothing for all sims):

      GenerateReport("Calibration"; sim_ids = ids, layout = :sequential,
                     snapshot_selection = :final)

LAYOUT (within each section, `layout` kwarg)
  :auto       — detect from the params that vary in the section:
                  0 varied → sequential, 1 varied → sorted 1D grid,
                  2 varied forming a full grid (GridVariation) → 2D tableau,
                  otherwise (scattered 2D / LHS, or 3+ params) → sequential grid.
  :sequential — ordered snapshots, in the order given (e.g. SMC posteriors).
  :grid       — force a 2D tableau on the first two varied params (sparse cells → N/A).

OUTPUT
  <output_dir>/index.html      the dashboard
  <output_dir>/snapshots/*.svg the copied snapshots (kept out of the top level)

Dependencies: PhysiCellModelManager (v0.3+), Dates, Printf
              BM metrics from ParameterOptimization.jl (included automatically)
=#

using PhysiCellModelManager
using Dates
using Printf

# BM morphology metrics (IC, area_frac_change, max_displacement, ...)
include(joinpath(@__DIR__, "ParameterOptimization.jl"))

# HTML template helpers (report_css/js, build_card_html, build_grid_html, build_2d_table_html)
include(joinpath(@__DIR__, "report_html.jl"))

# ══════════════════════════════════════════════════════════════════
#                    CONFIGURATION
# ══════════════════════════════════════════════════════════════════

Base.@kwdef struct ReportConfig
    output_dir::String
    sections::Vector{Pair{String,Any}} = Pair{String,Any}[]  # "Label" => ids|output, one section each
    sim_ids::Vector{Int} = Int[]                # used only when `sections` is empty (single section)
    layout::Symbol = :auto                      # :auto | :sequential | :grid
    min_cell_count::Int = 0                     # 0 = no filtering
    ignore_sims::Vector{Int} = Int[]
    param_aliases::Dict{String,String} = Dict{String,String}()
    snapshot_selection::Symbol = :by_cell_count # :final | :by_cell_count | :by_index
    target_cell_count::Union{Nothing,Int} = nothing
    snapshot_index::Int = -1                    # for :by_index mode
end

# Per-simulation record assembled for rendering.
const SimEntry = NamedTuple{(:dest_name, :cell_count, :param_vals, :bm_metrics),
                            Tuple{String, Int, Dict{String,String}, Dict{String,Float64}}}

# ══════════════════════════════════════════════════════════════════
#                    SNAPSHOT SELECTION
# ══════════════════════════════════════════════════════════════════

"""Construct the deterministic SVG path for a given simulation and snapshot index."""
function snapshot_svg_path(sim_id::Int, index)::String
    out_dir = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")
    index === :final && return joinpath(out_dir, "final.svg")
    return joinpath(out_dir, @sprintf("snapshot%08d.svg", index))
end

"""
    select_snapshot(sim_id, config) → (svg_path, cell_count)

Pick the snapshot for a simulation. Cell counts come from PCMM's population time
series rather than from parsing SVG circles.
"""
function select_snapshot(sim_id::Int, config::ReportConfig)::Tuple{String,Int}
    out_dir = joinpath(PhysiCellModelManager.dataDir(), "outputs", "simulations", "$sim_id", "output")
    !isdir(out_dir) && return ("", 0)

    final_count() = sum(values(try finalPopulationCount(sim_id) catch; Dict{String,Int}() end); init=0)

    if config.snapshot_selection == :final
        p = snapshot_svg_path(sim_id, :final)
        return (isfile(p) ? p : "", final_count())

    elseif config.snapshot_selection == :by_cell_count
        target = something(config.target_cell_count, 0)
        # Shared with the calibration path (ParameterOptimization.jl) so the
        # report and the fit metrics select the same frame.
        res = try
            output_index_at_cell_count(sim_id, target)
        catch e
            @warn "Could not load population time series for sim $sim_id" exception=e
            nothing
        end
        if res !== nothing
            idx, total = res
            p = snapshot_svg_path(sim_id, idx)   # 0-based output index
            isfile(p) && return (p, total)
        end
        p = snapshot_svg_path(sim_id, :final)  # fallback to final
        return (isfile(p) ? p : "", final_count())

    elseif config.snapshot_selection == :by_index
        svgs = sort(filter(f -> startswith(f, "snapshot") && endswith(f, ".svg"), readdir(out_dir)))
        isempty(svgs) && return ("", 0)
        idx = config.snapshot_index < 0 ? length(svgs) + config.snapshot_index + 1 : config.snapshot_index + 1
        idx = clamp(idx, 1, length(svgs))
        return (joinpath(out_dir, svgs[idx]), final_count())
    end
    return ("", 0)
end

# ══════════════════════════════════════════════════════════════════
#                    SMALL HELPERS
# ══════════════════════════════════════════════════════════════════

"""Display name for a parameter (apply alias, strip PCMM path prefixes)."""
function display_name(param::String, aliases::Dict{String,String})::String
    haskey(aliases, param) && return aliases[param]
    clean = replace(param, "user_parameters/" => "", "user_parameter/" => "")
    return clean
end

"""Format a numeric value for display (drop trailing zeros)."""
function format_val(v)::String
    v === missing && return "—"
    if v isa AbstractFloat
        (v == floor(v) && abs(v) < 1e10) && return string(Int(v))
        return string(round(v; sigdigits=4))
    end
    return string(v)
end

# Parse a formatted value for numeric sorting (non-numbers sort first).
numval(v::AbstractString) = (x = tryparse(Float64, v); x === nothing ? 0.0 : x)

# Make a string safe to use inside a filename.
sanitize(s::AbstractString) = replace(s, r"[^\w.\-]+" => "_")

"""Sorted, unique non-empty values of param `p` across `sids`."""
function param_values(sids, sim_data, p)
    vs = [get(sim_data[s].param_vals, p, "") for s in sids]
    sort(unique(filter(!isempty, vs)); by=numval)
end

"""True when the sims tile a complete `nx × ny` grid (one sim per cell, no gaps)."""
is_full_grid(xvals, yvals, present::Set, n) =
    n == length(xvals) * length(yvals) && length(present) == n

"""Resolve a `sections` value (Vector{Int} or a run/trial output) to a Vector{Int}."""
to_sim_ids(v)::Vector{Int} = v isa AbstractVector{<:Integer} ? collect(Int, v) : collect(Int, simulationIDs(v))

# ══════════════════════════════════════════════════════════════════
#                    SECTION RENDERING
# ══════════════════════════════════════════════════════════════════

"""
    render_section(label, idx, sids, sim_data, varied, layout, aliases) → (nav, body)

Render one section as HTML, choosing 2D tableau vs. grid-of-cards from `layout`
and the params that vary within the section (`varied`).
"""
function render_section(label, idx, sids, sim_data, varied, layout, aliases)
    nav, body = "", ""
    if !isempty(label)
        nav  *= """<li class="section-header">$(label)</li>\n"""
        body *= """<div class="exp-section" id="exp-$(idx)">$(label)</div>\n"""
    end

    # --- Decide whether to draw a 2D tableau, and on which two params ---
    table_params = String[]
    if layout == :grid && length(varied) >= 2
        table_params = sort(varied)[1:2]
    elseif layout == :auto && length(varied) == 2
        py, px = sort(varied)
        present = Set((get(sim_data[s].param_vals, py, ""), get(sim_data[s].param_vals, px, "")) for s in sids)
        if is_full_grid(param_values(sids, sim_data, px), param_values(sids, sim_data, py), present, length(sids))
            table_params = [py, px]
        end
    end

    anchor = "s-$(idx)"
    if !isempty(table_params)
        # ═════ 2D TABLE ═════
        py, px = table_params
        dy, dx = display_name(py, aliases), display_name(px, aliases)
        nav  *= """<li><a href="#$(anchor)">$(dy) × $(dx)</a></li>\n"""
        body *= """<h2 id="$(anchor)">$(dy) × $(dx)</h2>\n"""

        all_y, all_x = param_values(sids, sim_data, py), param_values(sids, sim_data, px)
        matrix = Dict{Tuple{String,String},Tuple{String,Int,Dict{String,Float64}}}()
        for s in sids
            d = sim_data[s]
            yv, xv = get(d.param_vals, py, ""), get(d.param_vals, px, "")
            (yv == "" || xv == "") && continue
            matrix[(yv, xv)] = (d.dest_name, d.cell_count, d.bm_metrics)
        end
        body *= build_2d_table_html(dy, dx, all_y, all_x, matrix)
    else
        # ═════ GRID OF CARDS (1D / sequential / 3+) ═════
        heading = isempty(varied) ? (isempty(label) ? "Snapshots" : label) :
                  join([display_name(p, aliases) for p in sort(varied)], " + ") *
                  (length(varied) == 1 ? " Variations" : "")
        nav  *= """<li><a href="#$(anchor)">$(heading)</a></li>\n"""
        body *= """<h2 id="$(anchor)">$(heading)</h2>\n"""

        # :auto with varied params sorts by the primary one; :sequential keeps input order.
        ordered = (layout != :sequential && !isempty(varied)) ?
                  sort(sids; by=s -> numval(get(sim_data[s].param_vals, sort(varied)[1], ""))) : sids

        cards = String[]
        for s in ordered
            d = sim_data[s]
            lbl = isempty(varied) ? "Sim $(s)" :
                  join(["$(display_name(p, aliases))=$(get(d.param_vals, p, ""))" for p in sort(varied)], ", ")
            push!(cards, build_card_html(d.dest_name, lbl, d.cell_count; metrics=d.bm_metrics))
        end
        body *= build_grid_html(cards)
    end
    return nav, body
end

# ══════════════════════════════════════════════════════════════════
#                    MAIN FUNCTION
# ══════════════════════════════════════════════════════════════════

function GenerateReport(output_dir::String; kwargs...)
    config = ReportConfig(; output_dir=output_dir, kwargs...)

    snap_dir = joinpath(config.output_dir, "snapshots")
    mkpath(snap_dir)
    println(" GenerateReport starting...")
    println("   Output: $(config.output_dir)")

    # --- Resolve ordered sections: (label, sim_ids) ---
    drop_ignored(ids) = filter(id -> !(id in config.ignore_sims), ids)
    if !isempty(config.sections)
        raw_sections = [(label, drop_ignored(to_sim_ids(v))) for (label, v) in config.sections]
    else
        ids = isempty(config.sim_ids) ? sort(collect(simulationIDs())) : collect(config.sim_ids)
        raw_sections = [("", drop_ignored(ids))]
    end
    all_ids = unique(reduce(vcat, (ids for (_, ids) in raw_sections); init=Int[]))
    println("   Simulations: $(all_ids)")

    # --- Constant parameters (shown once at the top of the report) ---
    df_all = simulationsTable(all_ids; remove_constants=false)
    varied_all = filter(!=("SimID"), names(simulationsTable(all_ids)))
    display_constant_cols = filter(setdiff(names(df_all), ["SimID"; varied_all])) do col
        val = df_all[1, col]
        val !== missing && string(val) != ""
    end

    # --- Build per-section render data ---
    section_data = Tuple{String, Vector{Int}, Dict{Int,SimEntry}, Vector{String}}[]
    filtered_count = 0

    for (label, sids) in raw_sections
        isempty(sids) && continue

        # Params that vary within this section (needs >1 sim for PCMM to drop constants).
        varied = length(sids) > 1 ? filter(!=("SimID"), names(simulationsTable(sids))) : String[]
        df = simulationsTable(sids; remove_constants=false)

        sim_data = Dict{Int,SimEntry}()
        kept = Int[]
        for sid in sids
            svg_path, cc = select_snapshot(sid, config)
            svg_path == "" && continue
            if config.min_cell_count > 0 && cc < config.min_cell_count
                println("   Sim $(sid): FILTERED (cells=$(cc) < $(config.min_cell_count))")
                filtered_count += 1
                continue
            end

            # This sim's varied-parameter values.
            row = df[df.SimID .== sid, :]
            param_vals = Dict{String,String}()
            if size(row, 1) > 0
                for col in varied
                    param_vals[col] = format_val(row[1, col])
                end
            end

            # Unique destination filename (label + varied values) under snapshots/.
            prefix = isempty(label) ? "" : sanitize(label) * "_"
            stem = isempty(param_vals) ? "sim$(sid)" :
                   join(["$(display_name(p, config.param_aliases))_$(v)" for (p, v) in sort(collect(param_vals))], "_")
            dest_name = joinpath("snapshots", prefix * stem * ".svg")
            cp(svg_path, joinpath(config.output_dir, dest_name); force=true)

            bm = try evaluate_simulation(sid) catch; Dict{String,Float64}() end
            sim_data[sid] = (dest_name=dest_name, cell_count=cc, param_vals=param_vals, bm_metrics=bm)
            push!(kept, sid)

            ic = haskey(bm, "IC") ? @sprintf(" IC=%.3f", bm["IC"]) : ""
            pstr = join(["$(display_name(p, config.param_aliases))=$(v)" for (p, v) in sort(collect(param_vals))], ", ")
            println("   [$(isempty(label) ? "—" : label)] Sim $(sid): cells=$(cc), $(pstr)$(ic)")
        end
        isempty(kept) && continue
        push!(section_data, (label, kept, sim_data, varied))
    end

    println("   Sections: $(length(section_data)) | Copied: $(sum(length(s[2]) for s in section_data; init=0)) | Filtered: $(filtered_count)")

    # --- Assemble HTML ---
    timestamp = Dates.format(now(), "U d, yyyy \\a\\t HH:MM")
    const_items = join(
        ["<li><strong>$(display_name(c, config.param_aliases)):</strong> $(format_val(df_all[1, c]))</li>"
         for c in display_constant_cols], "\n")

    nav = """<div class="sidebar">
<button class="close-btn" onclick="toggleSidebar()" title="Close">&times;</button>
<h3>Navigation</h3>
<ul>\n"""
    body = """<div class="main">
<h1>Simulation Snapshot Comparisons</h1>
<div class="ts">Generated on $(timestamp)</div>
<div class="bb"><h4>Constant Parameters</h4><ul>$(const_items)</ul></div>\n"""

    if config.min_cell_count > 0 || !isempty(config.ignore_sims) || config.snapshot_selection == :by_cell_count
        body *= """<div class="filter-info"><h4>⚙ Active Filters</h4>\n"""
        config.min_cell_count > 0 && (body *= """<p>Min cell count: <strong>$(config.min_cell_count)</strong> ($(filtered_count) sim(s) excluded)</p>\n""")
        !isempty(config.ignore_sims) && (body *= """<p>Ignored sims: <strong>$(config.ignore_sims)</strong></p>\n""")
        body *= """<p>Snapshot selection: <strong>$(config.snapshot_selection)</strong>$(config.snapshot_selection == :by_cell_count ? " (target ≥ $(config.target_cell_count) cells)" : "")</p>\n</div>\n"""
    end

    for (idx, (label, sids, sim_data, varied)) in enumerate(section_data)
        sec_nav, sec_body = render_section(label, idx, sids, sim_data, varied, config.layout, config.param_aliases)
        nav *= sec_nav
        body *= sec_body
    end

    nav *= "</ul>\n</div>\n"
    body *= "</div>\n"

    full_html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><title>Model Comparisons Dashboard</title>
<style>$(report_css())</style></head>
<body>
<button class="open-btn" onclick="toggleSidebar()"><span style="font-size:18px">☰</span> Menu</button>
$(nav)$(body)$(report_js())
</body></html>"""

    report_path = joinpath(config.output_dir, "index.html")
    open(f -> write(f, full_html), report_path, "w")
    println(" Report generated: $(report_path)")
    return report_path
end
