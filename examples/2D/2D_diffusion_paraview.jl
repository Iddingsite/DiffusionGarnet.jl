using DiffusionGarnet
using DelimitedFiles
using Plots
using ProgressBars

cd(@__DIR__)
println("Number of threads: $(Threads.nthreads())")

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


time_save = collect(range(0, 1, length=51))u"Myr"  # define times at which to save

@unpack t_charact = domain2D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time

save_data_callback = PresetTimeCallback(time_save_ad, save_data_paraview)

path_save = "Grt_2D.h5"  # chose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)
sol = simulate(domain2D; callback=save_data_callback, path_save=path_save, save_everystep=false);
# 377.768768 seconds (16.12 M allocations: 15.863 GiB, 10.96% gc time, 5.41% compilation time)
