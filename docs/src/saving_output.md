# [Saving output as HDF5 files](@id saving_output)

When dealing with 2D or especially 3D data, it can be impractical or even impossible to keep every timestep in memory, as the RAM on a single machine can quickly become saturated. It may then be relevant to save certain timesteps of interest to disk for later post-processing. DiffusionGarnet has a built-in function for this purpose, using a callback function based on the [DiffEqCallbacks](https://docs.sciml.ai/DiffEqCallbacks/stable/) package to produce HDF5 files.

[HDF5](https://www.hdfgroup.org/solutions/hdf5/) is a data format designed to store and organise large amount of data and is the data format chosen for DiffusionGarnet.

As an example, we will use the data from the [Diffusion in 2D Cartesian coordinates on CPU](@ref 2D_diffusion_CPU) tutorial but this is also applicable to 1D Cartesian or spherical coordinates in the same manner.

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

We then need to decide at which timesteps we want to save our data, here at 0, 0.5 and 1 Myr:

```julia
time_save = [0, 0.5, 1]u"Myr"  # define times at which to save
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
3-element Vector{Float64}:
 0.0
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
path_save = "Grt_2D.h5"  # chose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)
sol = simulate(domain2D; callback=save_data_callback, path_save=path_save, save_everystep=false);
```

!!! note
    We use `save_everystep=false` in `simulation()` to prevent saving every timestep in the model as we save it to disk. That should speed up the computation and reduce the burden on the RAM.

This output:

```julia
Data saved at 0.0 Myr.
Data saved at 0.1 Myr.
Data saved at 0.2 Myr.
Data saved at 0.5 Myr.
Data saved at 1.0 Myr.
352.464431 seconds (17.47 M allocations: 1.789 GiB, 0.15% gc time, 4.01% compilation time)
```

which indicates the time at which we saved our data and the total run time of the simulation.

The output file produced can be read with the [official viewer](https://www.hdfgroup.org/downloads/hdfview/#download) or using [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) in Julia or other programming languages and external programs.

Using [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/), we can look at the structure of the HDF5 file produced:

```julia-repl
julia> using HDF5
julia> hf5 = h5open("Grt_2D.h5")
🗂️ HDF5.File: (read-only) Grt_2D.h5
└─ 📂 Diffusion_Grt
   ├─ 🏷️ CharacteristicDiffusionCoefficient
   ├─ 🏷️ CharacteristicLength
   ├─ 🏷️ CharacteristicTime
   ├─ 🏷️ Coordinates
   ├─ 🏷️ Dx(µm)
   ├─ 🏷️ Dy(µm)
   ├─ 🏷️ LengthX(µm)
   ├─ 🏷️ LengthY(µm)
   ├─ 🏷️ Nx
   ├─ 🏷️ Ny
   ├─ 🏷️ TotalTime(Myr)
   ├─ 📂 t0000
   │  ├─ 🏷️ Time(Myr)
   │  ├─ 📂 Ca
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Ca
   │  ├─ 📂 Fe
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Fe
   │  ├─ 📂 Mg
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Mg
   │  └─ 📂 Mn
   │     ├─ 🏷️ Center
   │     ├─ 🏷️ DataType
   │     └─ 🔢 Mn
   ├─ 📂 t0001
   │  ├─ 🏷️ CurrentDt(Myr)
   │  ├─ 🏷️ Time(Myr)
   │  ├─ 📂 Ca
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Ca
   │  ├─ 📂 Fe
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Fe
   │  ├─ 📂 Mg
   │  │  ├─ 🏷️ Center
   │  │  ├─ 🏷️ DataType
   │  │  └─ 🔢 Mg
   │  └─ 📂 Mn
   │     ├─ 🏷️ Center
   │     ├─ 🏷️ DataType
   │     └─ 🔢 Mn
   └─ 📂 t0002
      ├─ 🏷️ CurrentDt(Myr)
      ├─ 🏷️ Time(Myr)
      ├─ 📂 Ca
      │  ├─ 🏷️ Center
      │  ├─ 🏷️ DataType
      │  └─ 🔢 Ca
      ├─ 📂 Fe
      │  ├─ 🏷️ Center
      │  ├─ 🏷️ DataType
      │  └─ 🔢 Fe
      ├─ 📂 Mg
      │  ├─ 🏷️ Center
      │  ├─ 🏷️ DataType
      │  └─ 🔢 Mg
      └─ 📂 Mn
         ├─ 🏷️ Center
         ├─ 🏷️ DataType
         └─ 🔢 Mn
julia> close(hf5)  # close the HDF5 file
```

We can see that the HDF5 file holds some general information about the model, such as characteristic variables, the dimensions of the model and the three timesteps with the four major element compositions.