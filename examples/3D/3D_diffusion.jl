using DiffusionGarnet
using JLD2

cd(@__DIR__)
println("Number of threads: $(Threads.nthreads())")

# use JLD2 to load data
file = jldopen("3D_data.jld2", "r")
@unpack Mg0, Fe0, Mn0, Ca0, grt_boundary = file
close(file)

Lx = 9000.0u"µm"
Ly = 9000.0u"µm"
Lz = 9000.0u"µm"
tfinal = 1.0u"Myr"
T = 900u"°C"
P = 0.6u"GPa"

IC3D = InitialConditions3D(Mg0, Fe0, Mn0, Lx, Ly, Lz, tfinal; grt_boundary = grt_boundary)
domain3D = Domain(IC3D, T, P)

time_save = collect(range(0, 1, length=51))u"Myr"  # define times at which to save

@unpack t_charact = domain3D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time

save_data_callback = PresetTimeCallback(time_save_ad, save_data_paraview)

path_save = "Grt_3D.h5"  # chose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)

sol = simulate(domain3D; callback=save_data_callback, path_save=path_save, save_everystep=false);