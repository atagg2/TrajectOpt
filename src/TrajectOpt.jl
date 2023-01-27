module TrajectOpt

    using FLOWMath, SNOW, Plots, DifferentialEquations, 
    StaticArrays, DelimitedFiles, VortexLattice, CCBlade, LinearAlgebra
    DE = DifferentialEquations
    FM = FLOWMath
    SA = StaticArrays
    DF = DelimitedFiles
    VL = VortexLattice
    CC = CCBlade

    include("models.jl")
    include("dynamics.jl")
    include("forces.jl")
    include("optimize.jl")
    include("visualizations.jl")

    export Surface
    export Rotor
    export Environment 
    export Inertia
    export SimpleSurface
    export SimpleRotor
    export VLMSurface
    export CCBladeRotor
    export Parameters
    export Model

    export dynamics2D!
    export simulate
    export plot_simulation

    export OptimizationProblem
    export optimize
    export objective_constructor
    export trim_constraints
    export trajectory_constraint_constructor

    export forces_conventional_low_fidel
    export forces_conventional_mid_fidel
    export polar_constructor
    export solve_rotor

    export store_path
    export clear_paths
    export visualize_paths
    export truncate

end
