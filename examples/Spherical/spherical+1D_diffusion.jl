using DiffusionGarnet
using DelimitedFiles
using Plots
using Printf

cd(@__DIR__)

data = DelimitedFiles.readdlm("Data_Grt_Sph.txt", '\t', '\n', header=true)[1]

Mg0 = data[:, 4]
Fe0 = data[:, 2]
Mn0 = data[:, 3]
Ca0 = data[:, 5]
distance = data[:, 1]
Lr = Lx = (data[end,1] - data[1,1])u"µm"
tfinal = 15u"Myr"

ICSph = ICSphMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lr, tfinal)
IC1D = IC1DMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lx, tfinal)

T = 900u"°C"
P = 0.6u"GPa"

DomainSph = Domain(ICSph, T, P)
Domain1D = Domain(IC1D, T, P; bc_neumann = (true, false))

sol_sph = simulate(DomainSph)
sol_1D = simulate(Domain1D)

@unpack t_charact = DomainSph

anim = @animate for i = LinRange(0, sol_sph.t[end], 100)
    l = @layout [a ; b]

    p1 = plot(distance, Fe0, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = title = @sprintf("Total Time = %.2f Ma | T = %.0f °C | P = %.1f GPa", i*t_charact, T[1].val, P[1].val), legend=:outerbottomright, linecolor=1,xlabel = "Distance from the core (µm)")
    p1 = plot!(distance, sol_sph(i)[:,2], label="Fe Sph",linecolor=1, linewidth=1)
    p1 = plot!(distance, sol_1D(i)[:,2], label="Fe 1D",linecolor=5, linewidth=1)


    p2 = plot(distance, Mg0, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance from the core (µm)")
    p2 = plot!(distance, Mn0, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
    p2 = plot!(distance, Ca0, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
    p2 = plot!(distance, sol_sph(i)[:,1], label="Mg Sph",linecolor=2, linewidth=1)
    p2 = plot!(distance, sol_1D(i)[:,1], label="Mg 1D",linecolor=6, linewidth=1)

    p2 = plot!(distance, sol_sph(i)[:,3], label="Mn Sph", linecolor=3, linewidth=1)
    p2 = plot!(distance, sol_1D(i)[:,3], label="Mn 1D", linecolor=7, linewidth=1)

    p2 = plot!(distance, 1 .- sol_sph(i)[:,1] .- sol_sph(i)[:,2] .- sol_sph(i)[:,3], label="Ca Sph", linecolor=4, linewidth=1)
    p2 = plot!(distance, 1 .- sol_sph(i)[:,1] .- sol_sph(i)[:,2] .- sol_sph(i)[:,3], label="Ca 1D", linecolor=8, linewidth=1)

    plot(p1, p2, layout = l)
end every 1

println("Now, generating the gif...")
gif(anim, "Grt_Spherical+1D.gif", fps = 7)
println("...Done!")
