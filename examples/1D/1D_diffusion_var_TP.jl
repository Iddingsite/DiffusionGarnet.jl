using DiffusionGarnet
using DelimitedFiles
using Plots
using Printf

cd(@__DIR__)

const data = DelimitedFiles.readdlm("Data_Grt_1D.txt", '\t', '\n', header=true)[1]

const Mg0 = data[:, 4]
const Fe0 = data[:, 2]
const Mn0 = data[:, 3]
const Ca0 = data[:, 5]
const distance = data[:, 1]
const Lx = (data[end,1] - data[1,1])u"µm"
const tfinal = 15u"Myr"

IC1D = IC1DMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lx=Lx, tfinal=tfinal)

time_update = [0, 2, 4, 6, 8, 10, 12, 14]u"Myr"
T = [900, 850, 800, 750, 700, 650, 600, 550]u"°C"
P = [0.6, 0.5, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3]u"GPa"

domain1D = Domain(IC1D, T, P, time_update)

@unpack time_update_ad = domain1D

update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

sol = simulate(domain1D; callback=update_diffusion_coef_call);

# plotting
anim = @animate for i = LinRange(0, sol.t[end], 100)
    l = @layout [a ; b]

    p1 = plot(distance, Fe0, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
    p1 = plot!(distance, sol(i)[:,2], label="Fe",linecolor=1, linewidth=1)

    # plot title with pressure and temperature
    j = findfirst(x -> i < x, time_update_ad)
    if j !== nothing
        title!(p1, @sprintf("Total Time = %.2f Ma | T = %.0f °C | P = %.1f GPa", i, T[j].val, P[j].val))
    else
        title!(p1, @sprintf("Total Time = %.2f Ma | T = %.0f °C | P = %.1f GPa", i, T[end].val, P[end].val))
    end

    p2 = plot(distance, Mg0, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
    p2 = plot!(distance, Mn0, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
    p2 = plot!(distance, Ca0, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
    p2 = plot!(distance, sol(i)[:,1], label="Mg",linecolor=2, linewidth=1)

    p2 = plot!(distance, sol(i)[:,3], label="Mn", linecolor=3, linewidth=1)

    p2 = plot!(distance, 1 .- sol(i)[:,1] .- sol(i)[:,2] .- sol(i)[:,3], label="Ca", linecolor=4, linewidth=1)

    plot(p1, p2, layout = l)
end every 1

println("Now, generating the gif...")
gif(anim, "Grt_1D_var_TP.gif", fps = 7)
println("...Done!")
