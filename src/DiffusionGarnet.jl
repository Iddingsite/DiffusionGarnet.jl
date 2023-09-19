module DiffusionGarnet

using Reexport

@reexport using Plots
@reexport using OrdinaryDiffEq
@reexport using BenchmarkTools
@reexport using Symbolics
@reexport using Parameters
@reexport using Unitful
@reexport using DelimitedFiles
@reexport using LinearAlgebra
@reexport using DiffEqCallbacks
# @reexport using ParallelStencil
@reexport using Logging: global_logger
@reexport using TerminalLoggers: TerminalLogger

function __init__()
    # initialise global logger for OrdinaryDiffEq
    global_logger(TerminalLogger())
end

include("input/initialconditions.jl")
include("callbacks/update_diffusion_coef_TP.jl")
include("Discretisation/1D/semi_discretisation_1D.jl")
include("Discretisation/Spherical/semi_discretisation_spherical.jl")
include("simulate/simulate.jl")

export InitialConditions1D, InitialConditions2D, InitialConditions3D, InitialConditionsSpherical
export D_ini!, Domain
export semi_discretisation_diffusion_1D, semi_discretisation_diffusion_spherical
export simulate
export update_diffusion_coef

end # module DiffusionGarnet
