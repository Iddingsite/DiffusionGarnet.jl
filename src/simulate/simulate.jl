
function simulate(domain::Domain1D; callbacks=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_1D, u0, t, domain)

    if callbacks === nothing
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)
    else
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callbacks)
    end

    return sol
end

function simulate(domain::DomainSpherical; callbacks=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_spherical, u0, t, domain)

    if callbacks === nothing
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)
    else
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callbacks)
    end

    return sol
end