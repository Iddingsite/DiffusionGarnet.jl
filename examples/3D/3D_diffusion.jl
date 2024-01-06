using DiffusionGarnet
using NPZ
using Plots
using ProgressBars

cd(@__DIR__)
println("Number of threads: $(Threads.nthreads())")


Mg0 = npzread("3D_Data/CMg_LR.npy")
Fe0 = npzread("3D_Data/CFe_LR.npy")
Mn0 = npzread("3D_Data/CMn_LR.npy")
Ca0 = npzread("3D_Data/CCa_LR.npy")
grt_boundary = npzread("3D_Data/sphere_contour_LR.npy")


Lx = 9000.0u"µm"
Ly = 9000.0u"µm"
Lz = 9000.0u"µm"
tfinal = 1.0u"Myr"
T = 900u"°C"
P = 0.6u"GPa"

grt_boundary = convert(Array{Float64,3}, grt_boundary)

IC3D = InitialConditions3D(Mg0, Fe0, Mn0, Lx, Ly, Lz, tfinal; grt_boundary = grt_boundary)
domain3D = Domain(IC3D, T, P)

time_save = collect(range(0, 1, length=51))u"Myr"  # define times at which to save

@unpack t_charact = domain3D  # unpack characteristic time to nondimensionalise the time for the simulation
time_save_ad = ustrip.(u"Myr", time_save) ./ t_charact  # convert to Myr, remove units, and convert to nondimensional time

save_data_callback = PresetTimeCallback(time_save_ad, save_data_paraview)

path_save = "Grt_3D.h5"  # chose the name and the path of the HDF5 output file (make sure to add .h5 or .hdf5 at the end)

sol = simulate(domain3D; callback=save_data_callback, path_save=path_save, save_everystep=false);


@unpack tfinal_ad, t_charact = domain3D
distance = LinRange(0, ustrip(u"µm", Lx), size(Mg0,1))

println("Plotting...")
anim = @animate for i = tqdm(LinRange(0, tfinal_ad, 20))

    time = round(i*t_charact;digits=2)

    l = @layout [a b ; c d ]
    Ca = 1 .- sol(i)[:,:,1] .- sol(i)[:,:,2] .- sol(i)[:,:,3]
    replace!(Ca, 1=>0)

    p1 = heatmap(distance, distance, sol(i)[:,:,1], label="Mg", dpi=200, title="Mg", clim=(0, maximum(sol(0)[:,:,1])), ylabel= "Distance (µm)")
    p2 = heatmap(distance, distance, sol(i)[:,:,2], label="Fe", dpi=200, title="Fe", clim=(0.7, maximum(sol(0)[:,:,2])))
    p3 = heatmap(distance, distance, sol(i)[:,:,3], label="Mn", dpi=200, title="Mn", clim=(0, maximum(sol(0)[:,:,3])), xlabel= "Distance (µm)", ylabel= "Distance (µm)")
    p4 = heatmap(distance, distance, Ca, label="Ca", dpi=200, title="Ca", clim=(0, 0.1), xlabel= "Distance (µm)")

    plot(p1, p2, p3, p4, layout = l , plot_title="Total Time = $(time) Ma")
end every 1

println("...Done!")

println("Now, generating the gif...")
gif(anim, "Grt_2D_900.gif", fps = 3)
println("...Done!")