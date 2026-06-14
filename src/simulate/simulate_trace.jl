
"""
    simulate(domain::Domain1DTrace; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-6, kwargs...)

Solve the linear diffusion equation in 1D. Save all timesteps in the output solution type variable by default.

"""
function simulate(domain::Domain1DTrace; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-6, p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    (; tfinal_ad, u0) = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian_trace, u0, t, p)

    sol = @time solve(prob, solver; abstol=abstol, reltol=reltol, kwargs...)

    return sol
end

"""
    simulate(domain::DomainSphericalTrace; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-6, p_extra=nothing, kwargs...)

Solve the linear diffusion equation for trace elements in spherical coordinates. Save all timesteps in the output solution type variable by default.
"""
function simulate(domain::DomainSphericalTrace; path_save=nothing, solver=ROCK2(), abstol=1e-8, reltol=1e-6, p_extra=nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    (; tfinal_ad, u0) = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_spherical_trace, u0, t, p)

    sol = @time solve(prob, solver; abstol=abstol, reltol=reltol, kwargs...)

    return sol
end

"""
    simulate(domain::Domain2DTrace; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

Solve the linear diffusion equations for trace elements in 2D. Save only the first and last timestep in the output solution type variable by default.
"""
function simulate(domain::Domain2DTrace; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    (; tfinal_ad, u0) = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian_trace, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    return sol
end

"""
    simulate(domain::Domain3DTrace; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

Solve the linear diffusion equations in 3D. Save only the first and last timestep in the output solution type variable by default.
"""
function simulate(domain::Domain3DTrace; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

    p = (domain = domain, path_save = path_save, p_extra = p_extra)

    (; tfinal_ad, u0) = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_cartesian_trace, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    return sol
end
