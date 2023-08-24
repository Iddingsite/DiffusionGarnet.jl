
function simulate(domain::Domain1D; callbacks=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_dicretisation_diffusion_1D, u0, t, domain)

    @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true)

    return sol
end