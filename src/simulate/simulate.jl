
"""
    simulate(domain; path_save=nothing, kwargs...)

Solve the coupled major element diffusion equations for a given domain using finite differences for the discretisation in space and return a solution type variable.

The time discretisation is based on the ROCK2 method, a stabilized explicit method (Adbdulle and Medovikov, 2001 ; https://doi.org/10.1007/s002110100292) using OrdinaryDiffEq.jl.

The solution type variable is following the format of OrdinaryDiffEq.jl (see https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/), and can be used to plot the solution, and to extract the solution at a given time. As the system is nondimensionalised, the time of the solution is in nondimensional time.

`path_save` is an optional argument, which can be used to define the path of the HDF5 output file. Default is to nothing.

All other accepted arguments such as `callback` or `progress` are the same as those of the `solve` function from DifferentialEquations (see https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).
"""
function simulate end

"""
    simulate(domain::Domain1D; callback=nothing, path_save=nothing, progress=true, save_everystep=true)

Solve the coupled major element diffusion equations in 1D. Save all timesteps in the output solution type variable by default.

"""
function simulate(domain::Domain1D; path_save=nothing, solver=ROCK2(), kwargs...)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0, t_charact = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_1D, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    sol.t .*= t_charact

    return sol
end

"""
    simulate(domain::DomainSpherical; callback=nothing, path_save=nothing, progress=true, save_everystep=true)

Solve the coupled major element diffusion equations in spherical coordinates. Save all timesteps in the output solution type variable by default.
"""
function simulate(domain::DomainSpherical; path_save=nothing, solver=ROCK2(), kwargs...)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0, t_charact = p.domain

    t = (0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_spherical, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    sol.t .*= t_charact

    return sol
end

"""
    simulate(domain::Domain2D; callback=nothing, path_save=nothing, progress=true, save_everystep=false, save_start=true)

Solve the coupled major element diffusion equations in 2D. By default, save only the first and last timestep in the output solution type variable by default.
"""
function simulate(domain::Domain2D; path_save=nothing, solver=ROCK2(), kwargs...)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0, t_charact = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_2D, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    sol.t .*= t_charact

    return sol
end

"""
    simulate(domain::Domain3D; callback=nothing, path_save=nothing, progressbar=true, save_everystep=false, save_start=true)

Solve the coupled major element diffusion equations in 3D. Save only the first and last timestep in the output solution type variable by default.

"""
function simulate(domain::Domain3D; path_save=nothing, solver=ROCK2(), kwargs...)

    p = (domain = domain, path_save = path_save)

    @unpack tfinal_ad, u0, t_charact = p.domain

    t = (0.0, tfinal_ad)

    prob = ODEProblem(semi_discretisation_diffusion_3D, u0, t, p)

    sol = @time solve(prob, solver; kwargs...)

    sol.t .*= t_charact

    return sol
end
