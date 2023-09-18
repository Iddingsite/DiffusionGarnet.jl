
function simulate(domain::Domain1D; callbacks=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_1D, u0, t, domain)

    @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)

    return sol
end

function simulate(domain::DomainSpherical; callbacks=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_spherical, u0, t, domain)

    @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)

    return sol
end