using DiffusionGarnet

println(Threads.nthreads())

const USE_GPU = false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2);
else
    @init_parallel_stencil(Threads, Float64, 2);
end


cd(@__DIR__)

Mg0 = DelimitedFiles.readdlm("Xprp.txt", '\t', '\n', header=false)
Fe0 = DelimitedFiles.readdlm("Xalm.txt", '\t', '\n', header=false)
Mn0 = DelimitedFiles.readdlm("Xsps.txt", '\t', '\n', header=false)
Ca0 = DelimitedFiles.readdlm("Xgrs.txt", '\t', '\n', header=false)
grt_boundary = DelimitedFiles.readdlm("Contour_HR.txt", '\t', '\n', header=false)

Lx = 9000.0u"µm"
Ly = 9000.0u"µm"
tfinal = 1.0u"Myr"
T = 900u"°C"
P = 0.6u"GPa"
IC2D = InitialConditions2D(Mg0, Fe0, Mn0, Lx, Ly, tfinal; grt_boundary = grt_boundary)
domain2D = Domain(IC2D, T, P)

sol = simulate(domain2D)

@unpack tfinal_ad, t_charact = domain2D

distance = LinRange(0, ustrip(u"µm", Lx), size(Mg0,1))

using ProgressBars

println("Plotting...")
anim = @animate for i = tqdm(LinRange(0, tfinal_ad, 20))

    time = round(i*t_charact;digits=2)

    l = @layout [a b ; c d ]
    Ca_ini = 1 .- sol(0)[:,:,1] .- sol(0)[:,:,2] .- sol(0)[:,:,3]
    replace!(Ca_ini, 1=>0)
    Ca = 1 .- sol(i)[:,:,1] .- sol(i)[:,:,2] .- sol(i)[:,:,3]
    replace!(Ca, 1=>0)

    p1 = heatmap(distance, distance, sol(i)[:,:,1], label="Mg", dpi=200, title="Mg", clim=(0, maximum(sol(0)[:,:,1])), ylabel= "Distance (µm)")
    p2 = heatmap(distance, distance, sol(i)[:,:,2], label="Fe", dpi=200, title="Fe", clim=(0.7, maximum(sol(0)[:,:,2])))
    p3 = heatmap(distance, distance, sol(i)[:,:,3], label="Mn", dpi=200, title="Mn", clim=(0, maximum(sol(0)[:,:,3])), xlabel= "Distance (µm)", ylabel= "Distance (µm)")
    p4 = heatmap(distance, distance, Ca, label="Ca", dpi=200, title="Ca", clim=(0, 0.1), xlabel= "Distance (µm)")

    plot(p1, p2, p3, p4, layout = l , plot_title="Total Time = $(time) Ma, T=$(round(ustrip.(u"°C", T); digits=2)) °C")
end every 1

println("...Done!")

println("Now, generating the gif...")
gif(anim, "Grt_2D_900.gif", fps = 3)
println("...Done!")