using DiffusionGarnet

const data = DelimitedFiles.readdlm("./examples/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

const Mg_ini = data[:, 12]
const Fe_ini = data[:, 10]
const Mn_ini = data[:, 11]
const Ca_ini = data[:, 13]
const distance = data[:, 1]
const Lx = (data[end,1] - data[1,1])u"µm"
const tfinal = 0.2u"Myr"

l = @layout [a ; b]

p1 = plot(distance, Fe_ini, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Initial conditions", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")

p2 = plot(distance, Mg_ini, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
p2 = plot!(distance, Mn_ini, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
p2 = plot!(distance, Ca_ini, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)

display(plot(p1, p2, layout = l))

IC1D = InitialConditionsCartesian(Mg_ini, Fe_ini, Mn_ini, tfinal, Lx)