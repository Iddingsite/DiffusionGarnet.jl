module DiffusionGarnet

using Reexport: @reexport

@reexport using UnPack: @unpack
@reexport using Unitful
@reexport using DiffEqCallbacks: PresetTimeCallback
@reexport using OrdinaryDiffEqStabilizedRK: ROCK2, ROCK4, ESERK4, ESERK5
@reexport using GeoParams: SetChemicalDiffusion
@reexport using GeoParams.Garnet

using Logging: global_logger
using TerminalLoggers: TerminalLogger
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve
using ParallelStencil
using GeoParams: compute_D, SetChemicalDiffusion, AbstractChemicalDiffusion
using GeoParams.Garnet: Grt_Fe_Chakraborty1992, Grt_Mg_Chakraborty1992, Grt_Mn_Chakraborty1992
using GeoParams.Garnet: Grt_Fe_Carlson2006, Grt_Mg_Carlson2006, Grt_Mn_Carlson2006, Grt_Ca_Carlson2006
using GeoParams.Garnet: Grt_Fe_Chu2015, Grt_Mg_Chu2015, Grt_Mn_Chu2015, Grt_Ca_Chu2015
using Preferences
using HDF5: h5open, create_group, attributes, read_attribute
using DelimitedFiles
using TOML

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
    # thanks Pascal Aellig for this!
    if Base.PkgId(@__MODULE__) == Base.__toplevel__
        printstyled("""

            ▐▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▌
            ▐ ██████╗ ██╗███████╗███████╗██╗   ██╗███████╗██╗ ██████╗ ███╗   ██╗   ▌
            ▐ ██╔══██╗██║██╔════╝██╔════╝██║   ██║██╔════╝██║██╔═══██╗████╗  ██║   ▌
            ▐ ██║  ██║██║█████╗  █████╗  ██║   ██║███████╗██║██║   ██║██╔██╗ ██║   ▌
            ▐ ██║  ██║██║██╔══╝  ██╔══╝  ██║   ██║╚════██║██║██║   ██║██║╚██╗██║   ▌
            ▐ ██████╔╝██║██║     ██║     ╚██████╔╝███████║██║╚██████╔╝██║ ╚████║   ▌
            ▐ ╚═════╝ ╚═╝╚═╝     ╚═╝      ╚═════╝ ╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝   ▌
            ▐  ██████╗  █████╗ ██████╗ ███╗   ██╗███████╗████████╗     ██╗██╗      ▌
            ▐ ██╔════╝ ██╔══██╗██╔══██╗████╗  ██║██╔════╝╚══██╔══╝     ██║██║      ▌
            ▐ ██║  ███╗███████║██████╔╝██╔██╗ ██║█████╗     ██║        ██║██║      ▌
            ▐ ██║   ██║██╔══██║██╔══██╗██║╚██╗██║██╔══╝     ██║   ██   ██║██║      ▌
            ▐ ╚██████╔╝██║  ██║██║  ██║██║ ╚████║███████╗   ██║██╗╚█████╔╝███████╗ ▌
            ▐  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝╚═╝ ╚════╝ ╚══════╝ ▌
            ▐▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▌


        Version: $(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
        """, bold=true, color=:default)
    end

    # initialise global logger for OrdinaryDiffEq
    global_logger(TerminalLogger())

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

include("diffcoef/diffcoef.jl")
include("input/initialconditions.jl")
include("callbacks/update_diffusion_coef_TP.jl")
include("callbacks/output.jl")
include("callbacks/output_paraview.jl")
include("discretisation/1D/semi_discretisation_1D_major.jl")
# include("discretisation/1D/semi_discretisation_1D_trace.jl")
include("discretisation/2D/semi_discretisation_2D_major.jl")
include("discretisation/3D/semi_discretisation_3D_major.jl")
include("discretisation/spherical/semi_discretisation_spherical_major.jl")
# include("discretisation/spherical/semi_discretisation_spherical_trace.jl")
include("simulate/simulate.jl")

export IC1DMajor, IC2DMajor, IC3DMajor, ICSphMajor#, IC1DTrace
export Domain
export simulate
export update_diffusion_coef
export save_data, save_data_paraview

end
