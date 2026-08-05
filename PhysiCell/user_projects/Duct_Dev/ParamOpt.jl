#=
================================================================================
PCMM SYNTAX REFERENCE & EXAMPLES
================================================================================

ParamOpt.jl is the first script in the piepline for configuring our Ductal deformation PhysiCell simulations

PARAMETER VARIATIONS
   - DiscreteVariation (Fixed parameter value):
     DiscreteVariation(configPath("group", "param"), value)
     Example: DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1)

   - DistributedVariation (Range for parameter sweeps):
     DistributedVariation(configPath("group", "param"), Distribution)
     Example: DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.2))

LHS PARAMETER SWEEPS
   Runs simulations sampled across the parameter space.
   
   Syntax:
     run(inputs, variations...; sampler=LatinHypercube(N), n_replicates=R)
   Example:
     out = run(inputs, dv_max_time, 
         DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.2));
         sampler = LatinHypercube(25), n_replicates = 1
     )

COMPARATIVE REPORTS
   Generates interactive HTML dashboards comparing simulation metrics.
   
   Syntax:
     GenerateReport("Report_Name"; sim_ids=Vector{Int}, kwargs...)
     * snapshot_selection options: :final, :by_cell_count, :by_index
   Example:
     GenerateReport("Report_1_Strain_Force";
         sim_ids = simulationIDs(out),
         snapshot_selection = :by_cell_count,
         target_cell_count = 200
     )

================================================================================
=#

using PhysiCellModelManager
using Distributions # Required for Uniform() sampling
include("GenerateReport.jl")

setNumberOfParallelSims(10)

# Clear all previous simulations to ensure a clean slate for the reports
println("Clearing previous simulations database...")
resetDatabase(; force_reset=true, force_continue=true)

############ Setup ############

config_folder = custom_code_folder = rulesets_collection_folder = "DuctDev_ParamOpt" 
inputs = InputFolders(config_folder, custom_code_folder;
                      rulesets_collection=rulesets_collection_folder)

# Set max time to 5 days 
dv_max_time = DiscreteVariation(configPath("max_time"), 7200)
force_recompile = false

############ EXPERIMENT 1: Strain Force Sweep ############
println("Running Experiment 1: Strain Force (Linear & Exponential)...")

# 1A: Linear Strain Sweep
ref_strain_lin = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn ON Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 1)      # Use Linear
)
out_strain_lin = run(LHSVariation(30), ref_strain_lin,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.0005, 0.01));
    n_replicates = 1,
    force_recompile = force_recompile
)

# 1B: Exponential Strain Sweep
ref_strain_exp = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 1), # Turn ON Strain
    DiscreteVariation(configPath("user_parameters", "is_strain_lin"), 0)      # Use Exponential
)
out_strain_exp = run(LHSVariation(30), ref_strain_exp,
    DistributedVariation(configPath("user_parameters", "seg_lin"), Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters", "seg_exp"), Uniform(0.01, 0.05));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 2: Restoring Force Sweep ############
println("Running Experiment 2: Restoring Force (Linear & Exponential)...")

# 2A: Linear Restoring Force Sweep
ref_restore_lin = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 0), # Isolate restoring force
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 1)     # Use Linear
)
out_restore_lin = run(LHSVariation(30), ref_restore_lin,
    DistributedVariation(configPath("user_parameters", "home_lin"), Uniform(0.0005, 0.01));
    n_replicates = 1,
    force_recompile = force_recompile
)

# 2B: Exponential Restoring Force Sweep
ref_restore_exp = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "Segment_Elasticity"), 0), # Isolate restoring force
    DiscreteVariation(configPath("user_parameters", "is_restore_lin"), 0)     # Use Exponential
)
out_restore_exp = run(LHSVariation(40), ref_restore_exp,
    DistributedVariation(configPath("user_parameters", "home_lin"), Uniform(0.01, 0.1)),
    DistributedVariation(configPath("user_parameters", "home_exp"), Uniform(0.01, 0.05));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ EXPERIMENT 3: Lumenal Pressure Sweep ############
println("Running Experiment 3: Lumenal Pressure...")

ref_lumenal = createTrial(inputs, dv_max_time,
    DiscreteVariation(configPath("user_parameters", "is_lumenal_pressure"), 1) # Turn ON pressure
)
out_lumenal = run(LHSVariation(30), ref_lumenal,
    DistributedVariation(configPath("user_parameters", "lumenal_pressure_strength"), Uniform(0.0005, 0.01));
    n_replicates = 1,
    force_recompile = force_recompile
)

############ Generate Reports ############
println("Simulations complete. Generating specific reports...")

# Each run() becomes its own section (Linear vs Exponential).
# Report 1: Strain Force
GenerateReport("Reports/Report_1_Strain_Force";
    sections = ["Linear Strain"      => out_strain_lin,
                "Exponential Strain" => out_strain_exp],
    snapshot_selection = :by_cell_count,
    target_cell_count = 200
)

# Report 2: Restoring Force
GenerateReport("Reports/Report_2_Restoring_Force";
    sections = ["Linear Restoring"      => out_restore_lin,
                "Exponential Restoring" => out_restore_exp],
    snapshot_selection = :by_cell_count,
    target_cell_count = 200
)

# Report 3: Lumenal Pressure (single sweep → single section)
GenerateReport("Reports/Report_3_Lumenal_Pressure";
    sim_ids = simulationIDs(out_lumenal),
    snapshot_selection = :by_cell_count,
    target_cell_count = 200
)

println("All parameter sweeps and reports successfully generated!")