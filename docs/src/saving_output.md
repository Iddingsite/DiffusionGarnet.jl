# [Saving output as HDF5 files](@id saving_output)

When dealing with 2D or especially 3D data, it can be impractical or even impossible to keep every timestep in memory, as the RAM on a single machine can quickly become saturated. It may then be relevant to save certain timesteps of interest to disk for later post-processing. DiffusionGarnet has a built-in function for this purpose, using a callback function based on the [DiffEqCallbacks](https://docs.sciml.ai/DiffEqCallbacks/stable/) package to produce HDF5 files.

[HDF5](https://www.hdfgroup.org/solutions/hdf5/) is a data format designed to store and organise large amount of data. 

As an example, we will use the data from the [Diffusion in 2D Cartesian coordinates on CPU](@ref 2D_diffusion_CPU) tutorial.

!!! tip
    Make sure to start Julia with multiple threads of execution to speed up the simulation.

```julia
using DiffusionGarnet
using DelimitedFiles
using Plots

Mg0 = DelimitedFiles.readdlm("Xprp.txt", '\t', '\n', header=false)
Fe0 = DelimitedFiles.readdlm("Xalm.txt", '\t', '\n', header=false)
Mn0 = DelimitedFiles.readdlm("Xsps.txt", '\t', '\n', header=false)
Ca0 = DelimitedFiles.readdlm("Xgrs.txt", '\t', '\n', header=false)
grt_boundary = DelimitedFiles.readdlm("contour_Grt.txt", '\t', '\n', header=false)

Lx = 9000.0u"µm"
Ly = 9000.0u"µm"
tfinal = 1.0u"Myr"
T = 900u"°C"
P = 0.6u"GPa"

IC2D = InitialConditions2D(Mg0, Fe0, Mn0, Lx, Ly, tfinal; grt_boundary = grt_boundary)
domain2D = Domain(IC2D, T, P)
```

We then need to decide at which timesteps we want to save our data, here at 0, 0.1, 0.2, 0.5 and 1 Myr:

```julia
time_save = [0, 0.1, 0.2, 0.5, 1]u"Myr"  # define times at which to save
```

!!! note
    Make sure to always include 0 in your timesteps to save the initial conditions.

Two additional steps are then required, we need to convert this time to dimensionless time so that the solver can use it, and we need to define our callback function:

```julia
@unpack t_charact = domain2D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time
```

`time_save_ad` is the equivalent dimensionless time to time_save:

```julia
julia > time_save_ad
5-element Vector{Float64}:
 0.0
 0.00229522961746084
 0.00459045923492168
 0.011476148087304199
 0.022952296174608398
```

!!! warning
    The time provided to the callback function should always be dimensionless.

We can now define our callback function that will save the data at these timesteps, using [`PresetTimeCallback()`](https://docs.sciml.ai/DiffEqCallbacks/stable/timed_callbacks/#DiffEqCallbacks.PresetTimeCallback):

```julia
save_data_callback = PresetTimeCallback(time_save_ad, save_data)
```

`save_data` is the function that will be called at each timestep in our callback function. `save_data_callback` can then be passed to `simulate` after defining the path of the output file:

```julia
path_save = "Grt_2D.h5"  # chose the name and the path of the HDF5 output file
sol = simulate(domain2D; callback=save_data_callback, path_save=path_save, save_everystep=false)
```