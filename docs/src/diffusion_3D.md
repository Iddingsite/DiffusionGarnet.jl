# [Diffusion in 3D Cartesian coordinates on CPU](@id 3D_diffusion_CPU)

!!! note
    This tutorial is combining the knowledge acquired from the [2D modelling](@ref 2D_diffusion_CPU) and the [Saving output as HDF5 files](@ref saving_output) tutorials. Make sure to understand them before starting this one.

The 3D geometry of garnets can be obtained from computed tomography (CT) scans or other similar techniques. Combined with an initial composition, the full grain can be modelled with realistic geometry.

In this tutorial, we will use fake geometry data and fake composition data with a low resolution and to get a reasonable runtime. The acquisition and processing of original data is beyond the scope of this tutorial and this package, and is left to the user. The goal of this tutorial is only to illustrate the use of DiffusionGarnet.jl in 3D Cartesian coordinates and show how to visualise the results.

To do so, we will use a callback function to save the results of the simulation to disk at regular intervals to be able to visualize the results using the software [Paraview](https://www.paraview.org/).

As mentioned in the tutorial for [2D modelling](@ref 2D_diffusion_CPU) DiffusionGarnet internally uses the package [ParallelStencil.jl](https://github.com/omlins/ParallelStencil.jl). Make sure to start with multiple threads to get the most out of this approach if you run the model on CPU.

For this tutorial, we will use example data from the [3D examples section](https://github.com/Iddingsite/DiffusionGarnet.jl/tree/main/examples/3D), which contains four composition matrices, and the contour of a fake spherical grain. The total resolution is of 128x128x128 voxels.

First, we load the data, which should be in the same folder as your current session:

```julia
using DiffusionGarnet  # this can take a while
using JLD2
using Plots

println("Number of threads: $(Threads.nthreads())")
# should ideally output more than 1 thread

# use JLD2 to load data
file = jldopen("3D_data.jld2", "r")
@unpack Mg0, Fe0, Mn0, Ca0, grt_boundary = file
close(file)
```
We will use the same conditions as in the [2D modelling tutorial](@ref 2D_diffusion_CPU):

```julia
# define total length in x and y
Lx = 9000.0u"µm"
Ly = 9000.0u"µm"
Lz = 9000.0u"µm"
# define total time for the model
tfinal = 1.0u"Myr"
# define the pressure and temperature conditions
T = 900u"°C"
P = 0.6u"GPa"

IC3D = InitialConditions3D(Mg0, Fe0, Mn0, Lx, Ly, Lz, tfinal; grt_boundary = grt_boundary)
domain3D = Domain(IC3D, T, P)
```

Now, let's define a callback function to save the results of the simulation in a HDF5 file, similar to the tutorial [Saving output as HDF5 files](@ref saving_output):

```julia
@unpack t_charact = domain3D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time

save_data_callback = PresetTimeCallback(time_save_ad, save_data_paraview)

path_save = "Grt_3D.h5"  # choose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)
```

!!! note
    We used here `save_data_paraview()` instead of `save_data()` to save the data in a format that can be read by Paraview.

We can now use the function `simulate()` to solve our system:

```julia
# solve the problem using DifferentialEquations.jl
sol = simulate(Domain3D; callback=save_data_callback, path_save=path_save, save_everystep=false)
```

!!! note
    The calculation can take some time, about 4 minutes on my machine with 16 threads.

Two files should have been created in the same folder as your current session: `Grt_3D.h5` and `Grt_3D.xdmf`. The first one is the HDF5 file containing the results of the simulation, and the second one is an XDMF file describing the HDF5 file. This last file can be opened with visualisation software programs, such as [Paraview](https://www.paraview.org/).

!!! warning
    For any pieces of visualisation software, make sure you are opening the XDMF file and not the HDF5 file. For [Paraview](https://www.paraview.org/), select the `XDMF Reader` as the reader when you open your data. Only Paraview was tested with this package, but other software should work as well.




