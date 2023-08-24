module DiffusionGarnet

using Reexport

@reexport using Plots
@reexport using OrdinaryDiffEq, DiffEqCallbacks, LinearSolve
@reexport using BenchmarkTools
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

export InitialConditions1D, InitialConditions2D, InitialConditions3D
export D_ini!, Domain

end # module DiffusionGarnet
