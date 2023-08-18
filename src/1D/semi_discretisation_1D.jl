


@views av(A) = 0.5*(A[1:end-1].+A[2:end])
# macros to avoid array allocation
macro qx(D,C,ix, Δx_ad) esc(:((0.5*($D[$ix]+$D[$ix+1])) * ($C[$ix+1]-$C[$ix])/Δx_ad )) end
macro sum_D(CMg, CFe, CMn, D0, ix) esc(:($D0[1] * $CMg[$ix] + $D0[2] * $CFe[$ix] + $D0[3] * $CMn[$ix] +
                                       $D0[4] * (1 - $CMg[$ix] - $CFe[$ix] - $CMn[$ix]))) end










@parallel_indices (ix) function stencil_diffusion!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, Δx_ad)

    if ix<size(dtCMg,1) && ix > 1
        dtCMg[ix] = (@qx(DMgMg,CMg,ix,Δx_ad) - @qx(DMgMg,CMg,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DMgFe,CFe,ix,Δx_ad) - @qx(DMgFe,CFe,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DMgMn,CMn,ix,Δx_ad) - @qx(DMgMn,CMn,ix-1,Δx_ad)) / Δx_ad
        dtCFe[ix] = (@qx(DFeMg,CMg,ix,Δx_ad) - @qx(DFeMg,CMg,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DFeFe,CFe,ix,Δx_ad) - @qx(DFeFe,CFe,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DFeMn,CMn,ix,Δx_ad) - @qx(DFeMn,CMn,ix-1,Δx_ad)) / Δx_ad
        dtCMn[ix] = (@qx(DMnMg,CMg,ix,Δx_ad) - @qx(DMnMg,CMg,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DMnFe,CFe,ix,Δx_ad) - @qx(DMnFe,CFe,ix-1,Δx_ad)) / Δx_ad +
                    (@qx(DMnMn,CMn,ix,Δx_ad) - @qx(DMnMn,CMn,ix-1,Δx_ad)) / Δx_ad
    else
        dtCMg[ix] = 0
        dtCFe[ix] = 0
        dtCMn[ix] = 0
    end
    return
end


# time loop
function Euler_explicit!(CMg, CFe, CMn, dtCMg, dtCFe, dtCMn, p, t_ad)

    D, D0 = p
    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    Δt = Δx_ad^2/8
    nt = ceil(Int,t_ad[2] / Δt)
    nvis = 100

    for it = tqdm(1:nt)

        # update diffusive parameters
        @parallel Diffusion_para_1D!(CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, D0)

        # semi-discretization
        @parallel stencil_diffusion!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, Δx_ad)

        # Forward Euler
        CMg .+= dtCMg .* Δt
        CFe .+= dtCFe .* Δt
        CMn .+= dtCMn .* Δt

        if it%nvis == 0
            display(plot(distance,[CMg,CFe,CMn,1 .- CMg .- CFe .- CMn];xlims=(0,distance[end]), ylims=(0,1),
                          xlabel="distance", ylabel="Concentration",
                          title="time = $(round(it*Δt*t_charact,digits=1))"))
        end

    end
end

function semi_dicretisation_diffusion_Grt_DiffEq(du,u,p,t)

    D, D0 = p
    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    @parallel Diffusion_para_1D!(CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, D0)

    # semi-discretization
    @parallel stencil_diffusion!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, Δx_ad)
end


# Initialise 3d array of composition
u0 = @zeros((nx, 3))
CMg = @zeros(nx)
CFe = @zeros(nx)
CMn = @zeros(nx)

D0 = @zeros((4))
u0[:,1] .= Mg_ini
u0[:,2] .= Fe_ini
u0[:,3] .= Mn_ini
CMg .= Mg_ini
CFe .= Fe_ini
CMn .= Mn_ini

D0 = D_ini!(D0,T,P)

DMgMg = @zeros(nx)
DMgFe = @zeros(nx)
DMgMn = @zeros(nx)

DFeMg = @zeros(nx)
DFeFe = @zeros(nx)
DFeMn = @zeros(nx)

DMnMg = @zeros(nx)
DMnFe = @zeros(nx)
DMnMn = @zeros(nx)

dtCMg = @zeros(nx)
dtCFe = @zeros(nx)
dtCMn = @zeros(nx)

const D_charact = mean(D0)
const x_charact = L
const t_charact = x_charact^2 / D_charact
const Δx_ad = Δx / x_charact
# Initialise 3d array of interdiffusion coefficients D
@parallel Diffusion_para_1D!(CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, D0)
D = DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn


p = D, D0, sum
t_ad = [0.0, 5.0 / t_charact]  # time in Ma

# classic Forward Euler
@time Euler_explicit!(CMg, CFe, CMn, dtCMg, dtCFe, dtCMn, p, t_ad)
# 159.913832 seconds (225.35 M allocations: 36.822 GiB, 3.32% gc time)

CMg_Euler = copy(CMg)
CFe_Euler = copy(CFe)
CMn_Euler = copy(CMn)


# With DiffEq
prob = ODEProblem(semi_dicretisation_diffusion_Grt_DiffEq, u0, t_ad, p);
@time sol = solve(prob, TRBDF2(autodiff=false), progress = true, progress_steps = 1, save_start=true);
# 2.270926 seconds (1.28 M allocations: 153.180 MiB)

l = @layout [a ; b]

p1 = plot(distance, Fe_ini, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Timestep = $(t_ad[2] * t_charact) Ma", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
p1 = plot!(distance, sol[end][:,2], label="Fe DiffEq",linecolor=1)
p1 = plot!(distance, CFe_Euler, label="Fe Euler",linecolor=5, linestyle = :dot)

p2 = plot(distance, Mg_ini, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
p2 = plot!(distance, Mn_ini, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
p2 = plot!(distance, Ca_ini, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
p2 = plot!(distance, sol[end][:,1], label="Mg DiffEq",linecolor=2)
p2 = plot!(distance, CMg_Euler, label="Mg Euler",linecolor=6, linestyle = :dot)
p2 = plot!(distance, sol[end][:,3], label="Mn DiffEq", linecolor=3)
p2 = plot!(distance, CMn_Euler, label="Mn Euler",linecolor=7, linestyle = :dot)
p2 = plot!(distance, 1 .- sol[end][:,1] .- sol[end][:,2] .- sol[end][:,3], label="Ca DiffEq", linecolor=4)
p2 = plot!(distance, 1 .- CMg_Euler .- CFe_Euler .- CMn_Euler, label="Ca Euler", linecolor=8, linestyle = :dot)

plot(p1, p2, layout = l, dpi=200)
savefig("Grt_1D_comparaison.png")



anim = @animate for i = tqdm(LinRange(t_ad[1], t_ad[2], 200))
    l = @layout [a ; b]

    p1 = plot(distance, Fe_ini, label="Fe initial", linestyle = :dash, linewidth=1, dpi=200, title = "Timestep = $(i * t_charact) Ma", legend=:outerbottomright, linecolor=1,xlabel = "Distance (µm)")
    p1 = plot!(distance, sol(i)[:,2], label="Fe Correct",linecolor=1)

    p2 = plot(distance, Mg_ini, label="Mg initial", linestyle = :dash, linewidth=1, dpi=200,legend=:outerbottomright,linecolor=2,xlabel = "Distance (µm)")
    p2 = plot!(distance, Mn_ini, label="Mn initial", linestyle = :dash, linewidth=1, linecolor=3)
    p2 = plot!(distance, Ca_ini, label="Ca initial", linestyle = :dash, linewidth=1, linecolor=4)
    p2 = plot!(distance, sol(i)[:,1], label="Mg Correct",linecolor=2)
    p2 = plot!(distance, sol(i)[:,3], label="Mn Correct", linecolor=3)
    p2 = plot!(distance, 1 .- sol(i)[:,1] .- sol(i)[:,2] .- sol(i)[:,3], label="Ca Correct", linecolor=4)

    plot(p1, p2, layout = l)
end every 1


println("Now, generating the gif...")
gif(anim, "Grt_1D.gif", fps = 7)
println("...Done!")