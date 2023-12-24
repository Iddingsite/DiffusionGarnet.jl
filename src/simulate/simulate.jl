
"""
    simulate(domain; callback=nothing, progressbar=true)

Solve the coupled diffusion equation using finite differences for a given domain and return a solution type variable.

The time discretisation is based on the ROCK2 method, a stabilized explicit method (Adbdulle and Medovikov, 2001 ; https://doi.org/10.1007/s002110100292) using OrdinaryDiffEq.jl.

The solution type variable is following the format of OrdinaryDiffEq.jl (see https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/), and can be used to plot the solution, and to extract the solution at a given time. As the system is nondimensionalised, the time of the solution is in nondimensional time.

callback is an optional argument, which can be used to pass a callback function to the solver. It follows the format of DiffEqCallbacks.jl (see https://docs.sciml.ai/DiffEqCallbacks/stable/).

progressbar is an optional argument, which can be used to display a progressbar during the simulation. Default is to true.
"""
function simulate end

"""
    simulate(domain::Domain1D; callback=nothing, progressbar=true)

Solve the coupled diffusion equation in 1D. Save all timesteps in the output solution type variable.

"""
function simulate(domain::Domain1D; callback=nothing, path_save=nothing, progress=true)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0 = p.domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_1D, u0, t, p)

    sol = @time solve(prob, ROCK2(), progress=progress, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callback)

    return sol
end

"""
    simulate(domain::DomainSpherical; callback=nothing, progressbar=true)

Solve the coupled diffusion equation in spherical coordinates. Save all timesteps in the output solution type variable.

"""
function simulate(domain::DomainSpherical; callback=nothing, path_save=nothing, progress=true)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0 = p.domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_spherical, u0, t, p)

    sol = @time solve(prob, ROCK2(), progress=progress, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callback)

    return sol
end

"""
    simulate(domain::Domain2D; callback=nothing, progressbar=true)

Solve the coupled diffusion equation in 2D. Save only the first and last timestep in the output solution type variable.

"""
function simulate(domain::Domain2D; callback=nothing, path_save=nothing, progress=true, save_everystep=true)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0 = p.domain

    t = [0.0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_2D, u0, t, p)

    sol = @time solve(prob, ROCK2(), progress=progress, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, save_everystep = save_everystep, callback=callback)

    return sol
end
