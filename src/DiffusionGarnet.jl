module DiffusionGarnet

using Reexport

@reexport using Plots
@reexport using OrdinaryDiffEq, DiffEqCallbacks, LinearSolve
@reexport using BenchmarkTools
@reexport using Parameters
@reexport using Unitful
# @reexport using ParallelStencil
@reexport using Logging: global_logger
@reexport using TerminalLoggers: TerminalLogger

function __init__()
    # initialise global logger for DifferentialEquations
    global_logger(TerminalLogger())
end


include("input/geometry.jl")

export InitialConditions1D, InitialConditions2D, InitialConditions3D



end # module DiffusionGarnet
