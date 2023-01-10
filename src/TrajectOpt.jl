module TrajectOpt

    using FLOWMath, SNOW, Plots, DifferentialEquations, StaticArrays, DelimitedFiles, #=Snopt,=# VortexLattice
    DE = DifferentialEquations
    FM = FLOWMath
    SA = StaticArrays
    DF = DelimitedFiles
    # SN = Snopt
    VL = VortexLattice

    include("models.jl")
    include("dynamics.jl")
    include("optimize.jl")


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
    export forces_conventional_low_fidel
    export forces_conventional_mid_fidel

    export Model
    export LowFidel 
    export HighFidel
    export ConventionalLowFidel
    export ConventionalMidFidel
    export BiWingTailSitterLowFidel
    export polar_constructor
    export conventional_forces_constructor
    export biwing_tailsitter_forces_constructor
    export dynamics2D!
    export simulate
    export plot_simulation
    export optimize_trim
    export optimize_trajectory
    export optimize_trajectory_by_segments

end
