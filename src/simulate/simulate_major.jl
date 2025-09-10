

"""
    simulate(domain::Domain1DMajor; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-4, p_extra = nothing, kwargs...)

Solve the coupled major element diffusion equations in 1D. Save all timesteps in the output solution type variable by default.

"""
function simulate(domain::Domain1DMajor; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-4, p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    @unpack tfinal_ad, u0 = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian, u0, t, p)

    sol = @time solve(prob, solver; abstol=abstol, reltol=reltol, kwargs...)

    return sol
end

"""
    simulate(domain::DomainSphericalMajor; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-4, p_extra = nothing, kwargs...)

Solve the coupled major element diffusion equations in spherical coordinates. Save all timesteps in the output solution type variable by default.
"""
function simulate(domain::DomainSphericalMajor; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-4, p_extra=nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    @unpack tfinal_ad, u0 = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_spherical_major, u0, t, p)

    sol = @time solve(prob, solver; abstol=abstol, reltol=reltol, kwargs...)

    return sol
end

"""
    simulate(domain::Domain2DMajor; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

Solve the coupled major element diffusion equations in 2D. Save only the first and last timestep in the output solution type variable by default.
"""
function simulate(domain::Domain2DMajor; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    @unpack tfinal_ad, u0 = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    return sol
end

"""
    simulate(domain::Domain3DMajor; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

Solve the coupled major element diffusion equations in 3D. Save only the first and last timestep in the output solution type variable by default.

"""
function simulate(domain::Domain3DMajor; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    @unpack tfinal_ad, u0 = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    return sol
end
