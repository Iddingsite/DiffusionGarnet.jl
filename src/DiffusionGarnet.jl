module DiffusionGarnet

using Reexport

@reexport using Parameters
@reexport using Unitful
@reexport using DiffEqCallbacks
@reexport using Logging: global_logger
@reexport using TerminalLoggers: TerminalLogger
# @reexport using ParallelStencil

using OrdinaryDiffEq
using Symbolics
using ParallelStencil
using Preferences
using DelimitedFiles
using BenchmarkTools
using Plots

# initialise ParallelStencil (Thx AlbertDeMontserrat!)
function set_backend(new_backend::String)
    if !(
        new_backend ∈
        ("Threads_Float64_2D", "Threads_Float32_2D", "CUDA_Float64_2D", "CUDA_Float32_2D", "Threads_Float64_3D", "Threads_Float32_3D", "CUDA_Float64_3D", "CUDA_Float32_3D")
    )
        throw(ArgumentError("Invalid backend: \"$(new_backend)\""))
    end

    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("backend" => new_backend)
    @info("New backend set; restart your Julia session for this change to take effect!")
end

const backend = @load_preference("backend", "Threads_Float64_2D")

export backend, set_backend

function __init__()
    # initialise global logger for OrdinaryDiffEq
    global_logger(TerminalLogger())

    # @require ParallelStencil = "94395366-693c-11ea-3b26-d9b7aac5d958" include("Discretisation/2D/semi_discretisation_2D.jl")
end


let
    s = split(backend, "_")
    device = s[1]
    precision = s[2]
    dimension = parse(Int, s[3][1])
    @eval begin
        @init_parallel_stencil($(Symbol(device)), $(Symbol(precision)), $dimension)
    end
end

include("input/initialconditions.jl")
include("callbacks/update_diffusion_coef_TP.jl")
include("callbacks/output.jl")
include("Discretisation/1D/semi_discretisation_1D.jl")
include("Discretisation/2D/semi_discretisation_2D.jl")
include("Discretisation/Spherical/semi_discretisation_spherical.jl")
include("simulate/simulate.jl")

export InitialConditions1D, InitialConditions2D, InitialConditions3D, InitialConditionsSpherical
export D_ini!, Domain
export semi_discretisation_diffusion_1D, semi_discretisation_diffusion_spherical, semi_discretisation_diffusion_2D, Diffusion_coef_2D!, stencil_diffusion_2D!
export simulate
export update_diffusion_coef
export hdf5_initial_conditions, hdf5_timestep, save_data

end # module DiffusionGarnet