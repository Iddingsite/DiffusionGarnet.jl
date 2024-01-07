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
grt_boundary = DelimitedFiles.readdlm("grt_boundary.txt", '\t', Int, '\n', header=false)

Lx     = 9000.0u"µm"
Ly     = 9000.0u"µm"
tfinal = 1.0u"Myr"
T      = 900u"°C"
P      = 0.6u"GPa"

IC2D     = InitialConditions2D(Mg0, Fe0, Mn0, Lx, Ly, tfinal; grt_boundary = grt_boundary)
domain2D = Domain(IC2D, T, P)

sol = simulate(domain2D; save_everystep=true)
# 377.768768 seconds (16.12 M allocations: 15.863 GiB, 10.96% gc time, 5.41% compilation time)

@unpack tfinal_ad, t_charact = domain2D

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