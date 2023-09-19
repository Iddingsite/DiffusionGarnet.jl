using DiffusionGarnet

const data = DelimitedFiles.readdlm("./examples/1D/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

const Mg0 = data[:, 4]
const Fe0 = data[:, 2]
const Mn0 = data[:, 3]
const Ca0 = data[:, 5]
const distance = data[:, 1]
const Lx = (data[end,1] - data[1,1])u"µm"
const tfinal = 15u"Myr"

IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)

time_update = [0u"Myr", 2u"Myr", 4u"Myr", 6u"Myr", 8u"Myr", 10u"Myr", 12u"Myr", 14u"Myr"]
T = [900u"°C", 850u"°C", 800u"°C", 750u"°C", 700u"°C", 650u"°C", 600u"°C", 550u"°C"]
P = [0.6u"GPa", 0.5u"GPa", 0.4u"GPa", 0.3u"GPa", 3u"GPa", 3u"GPa", 3u"GPa", 3u"GPa"]
ustrip.(u"Myr", time_update)

domain1D = Domain(IC1D, T, P, time_update)

@unpack time_update_ad = domain1D

update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)


sol = simulate(domain1D; callbacks=update_diffusion_coef_call)


@unpack tfinal_ad, t_charact = domain1D

anim = @animate for i = LinRange(0, tfinal_ad, 100)
    l = @layout [a ; b]

    p1 = plot(distance, Fe0, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Timestep = $(round(((i)* t_charact);digits=2)) Ma", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
    p1 = plot!(distance, sol(i)[:,2], label="Fe",linecolor=1, linewidth=1)


    p2 = plot(distance, Mg0, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
    p2 = plot!(distance, Mn0, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
    p2 = plot!(distance, Ca0, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
    p2 = plot!(distance, sol(i)[:,1], label="Mg",linecolor=2, linewidth=1)

    p2 = plot!(distance, sol(i)[:,3], label="Mn", linecolor=3, linewidth=1)

    p2 = plot!(distance, 1 .- sol(i)[:,1] .- sol(i)[:,2] .- sol(i)[:,3], label="Ca", linecolor=4, linewidth=1)

    plot(p1, p2, layout = l)
end every 1

println("Now, generating the gif...")
gif(anim, "./examples/1D/Grt_1D_var_TP.gif", fps = 7)
println("...Done!")
