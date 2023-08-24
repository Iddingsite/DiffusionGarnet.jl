module DiffusionGarnet

using Reexport

@reexport using Plots
@reexport using OrdinaryDiffEq
@reexport using BenchmarkTools
@reexport using Symbolics
@reexport using Parameters
@reexport using Unitful
@reexport using DelimitedFiles
# @reexport using ParallelStencil
@reexport using Logging: global_logger
@reexport using TerminalLoggers: TerminalLogger

function __init__()
    # initialise global logger for OrdinaryDiffEq
    global_logger(TerminalLogger())
end

include("input/initialconditions.jl")
include("1D/semi_discretisation_1D.jl")
include("simulate/simulate.jl")

export InitialConditions1D, InitialConditions2D, InitialConditions3D
export D_ini!, Domain
export semi_dicretisation_diffusion_1D
export simulate

end # module DiffusionGarnet
