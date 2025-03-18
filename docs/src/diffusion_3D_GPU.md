# [Diffusion in 3D Cartesian coordinates on GPU](@id 3D_diffusion_GPU)

!!! note
    This tutorial is a follow-up to the [3D modelling](@ref 3D_diffusion_CPU) tutorial. Make sure to understand it before starting this one. Also, you need an NVIDIA GPU to run this tutorial.

To speed-up the simulation, it is possible to make use of GPU acceleration. This tutorial will show how to run the same simulation as in the [3D modelling](@ref 3D_diffusion_CPU) tutorial but on a GPU.

The first step consists as initialising DiffusionGarnet.jl to be compatible with GPUs using ParallelStencil.jl. To do so, we first need to load DiffusionGarnet.jl, modify its backend, and start again a new session.

To do so, you can run this command in your terminal in the folder containing your current project:

```
julia --project=. -e '
using DiffusionGarnet
set_backend("CUDA_Float32_3D")
exit()'
```

Going line by line, this command will:
1. Start Julia in the current project folder and activate the project environment.
2. Load DiffusionGarnet.jl.
3. Set the backend to `CUDA_Float32_3D`.
4. Exit Julia.

!!! note
    To switch back to multithread mode on CPU, replace `CUDA` by `Threads`.

Now, we can start a new Julia session and load DiffusionGarnet.jl. The simulation will now run on the GPU.

```julia
using DiffusionGarnet

println(backend)
````

This should now print: `CUDA_Float32_3D`.

The rest of the tutorial is the same as the [3D modelling](@ref 3D_diffusion_CPU) tutorial, except that the initial conditions and the domain need to be converted to `CuArray{Float32}` and not `Array{Float32}`:

```julia
using Downloads
using JLD2
using DiffusionGarnet
using CUDA

# define the current directory as the working directory
cd(@__DIR__)

# Define the Zenodo dataset URL, you can change the name to download other datasets in the Zenodo repository (https://zenodo.org/records/15045718)
# here, we will download the lowest resolution dataset (256³) to save time, for the model isolated matrix model (IMM). See publication for more details.
data_file = "256_cubed_IMM_compo.jld2"

# check if the file is already downloaded
if !isfile(data_file)
    # Define the Zenodo dataset URL, you can change the name to download other datasets in the Zenodo repository.
    zenodo_url = "https://zenodo.org/records/15045718/files/" * data_file * "?download=1"

    # Download the file in the same folder as this file (this can take a while if you connection is slow)
    Downloads.download(zenodo_url, data_file)
end

# use JLD2
file = jldopen(data_file, "r")
@unpack Mg0, Fe0, Mn0, Ca0, grt_boundary = file
close(file)

# define total length in x and y
Lx = 11422.61u"µm"
Ly = 11422.61u"µm"
Lz = 7623.57u"µm"
# define total time for the model
tfinal = 10.0u"Myr"
# define the pressure and temperature conditions
T = 700u"°C"
P = 0.8u"GPa"

# composition at the contact between garnet and matrix
Mg_border = 0.1152
Fe_border = 0.6012
Mn_border = 0.0435
Ca_border = 0.2401

# add this to fix the composition on the boundary
Mg0[grt_boundary .== 1] .= Mg_border
Fe0[grt_boundary .== 1] .= Fe_border
Mn0[grt_boundary .== 1] .= Mn_border
Ca0[grt_boundary .== 1] .= Ca_border

# convert to float32
Mg0 = convert(CuArray{Float32}, Mg0)
Fe0 = convert(CuArray{Float32}, Fe0)
Mn0 = convert(CuArray{Float32}, Mn0)
Ca0 = convert(CuArray{Float32}, Ca0)
grt_boundary = convert(CuArray{Float32}, grt_boundary)

IC3D = DiffusionGarnet.InitialConditions3D(Mg0, Fe0, Mn0, Lx, Ly, Lz, tfinal; grt_boundary=grt_boundary)
domain3D = Domain(IC3D, T, P)

# free memory, as 3D data is large
file = nothing
data = nothing
Mg0 = nothing
Fe0 = nothing
Mn0 = nothing
Ca0 = nothing
grt_boundary = nothing
IC3D = nothing

time_save_first = collect(range(0, 1, step=0.1))u"Myr"
time_save_second = collect(range(1.5, 10, step=0.5))u"Myr"
time_save = vcat(time_save_first,time_save_second)

@unpack t_charact = domain3D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time

# create the callback function
save_data_callback = PresetTimeCallback(time_save_ad, save_data_paraview, save_positions=(false,false))

path_save = "data_model_10_Ma_GPU.h5"  # chose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)

# run the simulation with ROCK2 solver
sol = simulate(domain3D; callback=save_data_callback, path_save=path_save, save_everystep=false,  save_start=false, progress=true, progress_steps=1, solver=ROCK2());
```