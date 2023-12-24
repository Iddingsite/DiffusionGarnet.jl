using DiffusionGarnet
using DelimitedFiles
using Plots

cd(@__DIR__)

data = DelimitedFiles.readdlm("Data_Grt_Sph.txt", '\t', '\n', header=true)[1]

Mg0 = data[:, 4]
Fe0 = data[:, 2]
Mn0 = data[:, 3]
Ca0 = data[:, 5]
distance = data[:, 1]
Lr = (data[end,1] - data[1,1])u"µm"
tfinal = 15u"Myr"

ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lr, tfinal)

T = 900u"°C"
P = 0.6u"GPa"

DomainSph = Domain(ICSph, T, P)

sol_sph = simulate(DomainSph)

@unpack tfinal_ad, t_charact = DomainSph

anim = @animate for i = LinRange(0, tfinal_ad, 100)
    l = @layout [a ; b]

    p1 = plot(distance, Fe0, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Total Time = $(round(((i)* t_charact);digits=2)) Ma", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
    p1 = plot!(distance, sol_sph(i)[:,2], label="Fe Sph",linecolor=1, linewidth=1)

    p2 = plot(distance, Mg0, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
    p2 = plot!(distance, Mn0, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
    p2 = plot!(distance, Ca0, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
    p2 = plot!(distance, sol_sph(i)[:,1], label="Mg Sph",linecolor=2, linewidth=1)

    p2 = plot!(distance, sol_sph(i)[:,3], label="Mn Sph", linecolor=3, linewidth=1)

    p2 = plot!(distance, 1 .- sol_sph(i)[:,1] .- sol_sph(i)[:,2] .- sol_sph(i)[:,3], label="Ca Sph", linecolor=4, linewidth=1)

    plot(p1, p2, layout = l)
end every 1

println("Now, generating the gif...")
gif(anim, "Grt_Spherical.gif", fps = 7)
println("...Done!")
