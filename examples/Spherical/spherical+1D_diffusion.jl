using DiffusionGarnet

const data = DelimitedFiles.readdlm("./examples/Spherical/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

const Mg0 = reverse(data[1:size(data,1)÷2, 4])
const Fe0 = reverse(data[1:size(data,1)÷2, 2])
const Mn0 = reverse(data[1:size(data,1)÷2, 3])
const Ca0 = reverse(data[1:size(data,1)÷2, 5])
const distance = data[1:size(data,1)÷2, 1]
const Lx = (data[size(data,1)÷2,1] - data[1,1])u"µm"
const tfinal = 15u"Myr"

ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lx, tfinal)

const T = 900u"°C"
const P = 0.6u"GPa"

DomainSph = Domain(ICSph, T, P)

IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)
Domain1D = Domain(IC1D, T, P; bc_neumann = (true, false))

sol_sph = simulate(DomainSph)
sol_1D = simulate(Domain1D)


@unpack tfinal_ad, t_charact = DomainSph

anim = @animate for i = LinRange(0, tfinal_ad, 100)
    l = @layout [a ; b]

    p1 = plot(distance, Fe0, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Total Time = $(round(((i)* t_charact);digits=2)) Ma", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
    p1 = plot!(distance, sol_sph(i)[:,2], label="Fe Sph",linecolor=1, linewidth=1)
    p1 = plot!(distance, sol_1D(i)[:,2], label="Fe 1D",linecolor=5, linewidth=1)


    p2 = plot(distance, Mg0, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
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
gif(anim, "./examples/Spherical/Grt_Spherical+1D.gif", fps = 7)
println("...Done!")
