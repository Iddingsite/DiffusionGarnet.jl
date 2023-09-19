
"""
    simulate(domain; callback=nothing)

Solve the coupled diffusion equation using finite differences for a given domain and return a solution type variable.

The time discretisation is based on the ROCK2 method, a stabilized explicit method (Adbdulle and Medovikov, 2001 ; https://doi.org/10.1007/s002110100292) using OrdinaryDiffEq.jl.

The solution type variable is following the format of OrdinaryDiffEq.jl (see https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/), and can be used to plot the solution, and to extract the solution at a given time. As the system is nondimensionalised, the time of the solution is in nondimensional time.

callbacks is an optional argument, which can be used to pass a callback function to the solver. It follows the format of DiffEqCallbacks.jl (see https://docs.sciml.ai/DiffEqCallbacks/stable/).
"""
function simulate end

"""
    simulate(domain::Domain1D; callback=nothing)

Solve the coupled diffusion equation in 1D. Save all timesteps in the output solution type variable.

"""
function simulate(domain::Domain1D; callback=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_1D, u0, t, domain)

    if callback === nothing
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)
    else
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callback)
    end

    return sol
end

"""
    simulate(domain::DomainSpherical; callback=nothing)

Solve the coupled diffusion equation in spherical coordinates. Save all timesteps in the output solution type variable.

"""
function simulate(domain::DomainSpherical; callback=nothing)

    @unpack tfinal_ad, u0 = domain

    t = [0, tfinal_ad]

    prob = ODEProblem(semi_discretisation_diffusion_spherical, u0, t, domain)

    if callback === nothing
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6)
    else
        @time sol = solve(prob, ROCK2(), progress=true, progress_steps=1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callback)
    end

    return sol
end