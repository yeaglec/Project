using PhysiCellModelManager
using Distributions
include("GenerateReport.jl") 

setNumberOfParallelSims(10)

println("=== Setting up ABC-SMC Calibration (Linear Mechanics) ===")
resetDatabase(; force_reset=true, force_continue=true)
println("Reset Database")

# Config
config_folder = "DuctDev_ParamOpt"
inputs = InputFolders(config_folder, config_folder; rulesets_collection=config_folder)

# Locking in the LINEAR mechanics
ref_model = createTrial(inputs,
    DiscreteVariation(configPath("max_time"), 11520), # 8 days
    
    # Strain Mechanics (ON, Linear)
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1),
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1),
    
    # Restoring Mechanics (ON, Linear)
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1),
    
    # Lumenal Pressure (ON)
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"), 1),

    # 4. Kernel Size (Static)
    DiscreteVariation(configPath("user_parameters", "membrane_force_smoothing_sigma"), 25.0)
)

# Evaluate descriptors at a fixed cell-count 
const TARGET_CELL_COUNT = 155
summary_at_milestone(monad_id) = bm_summary_statistic(monad_id; target_cell_count = TARGET_CELL_COUNT)

# Derived from your stable 7-day baseline
observed_target = Dict(
    "IC"                => 1.41,
    "area_frac_change"  => -0.095
)
# NEXT STEP (deferred ABC rework): widen the target to the normalized descriptors
# now emitted by evaluate_simulation, once per-metric weighting + a pilot-grounded
# epsilon are in place, e.g.:
#     "max_indent_depth" => ..., "indent_extent" => ..., "lobe_amp" => ...

############ 3. THE SEARCH SPACE (Priors) ############
# Tightly bounded based on your LHS empirical observations
parameters_to_tune = [
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters", "home_lin"), Uniform(0.0005, 0.005)),
    DistributedVariation(configPath("user_parameters", "lumenal_pressure_strength"), Uniform(0.005, 0.05))
]

############ 4. BUILD & RUN CALIBRATION ############
println("Building Calibration Problem...")
problem = CalibrationProblem(
    ref_model,
    parameters_to_tune,
    observed_target,
    summary_at_milestone,    # ← summary at the cell-count milestone
    bm_distance;
    n_replicates = 2
)

# Adaptive epsilon (next ε = median of accepted distances, floored at
# minimum_epsilon). The safety stops let the run end gracefully instead of
# grinding when ε can no longer be reduced.
println("Running ABC-SMC (Validating Run)...")
result = runABC(
    problem;
    population_size      = 50,
    max_nr_populations   = 6,
    minimum_epsilon      = 0.1,    # normalized-distance scale; tunable
    epsilon_quantile     = 0.5,
    min_acceptance_rate  = 0.02,   # stop if proposals stop getting accepted
    min_epsilon_decrease = 0.02,   # stop once ε plateaus
    description          = "Linear Mechanics Calibration (validating run)",
)

############ 5. EXTRACT & REPORT ############
df, weights = posterior(result)
println("\n═══ Final Posterior Parameter Estimates ═══")
println(df)
println("\n═══ Convergence (ε / acceptance_rate / ess per generation) ═══")
println(ConvergenceSummary(result))

# posterior() returns only parameter columns — get sim IDs from the calibration's monads.
# calibrationMonadIDs is not re-exported, so reach it through the ModelManager submodule.
monad_ids = PhysiCellModelManager.ModelManager.calibrationMonadIDs(result.calibration)
calibration_sim_ids = sort!(unique!(reduce(vcat, simulationIDs.(Monad.(monad_ids)); init = Int[])))

println("\n═══ Generating Report for Best Fits ═══")
GenerateReport("Reports/Calibration_Results_Linear_Mechanics";
    sim_ids = calibration_sim_ids,
    layout = :sequential,          # one snapshot per posterior parameter set, in order
    snapshot_selection = :final
)

