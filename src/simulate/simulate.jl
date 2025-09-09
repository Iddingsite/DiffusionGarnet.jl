"""
    simulate(domain; path_save=nothing, solver=ROCK2(), p_extra = nothing, kwargs...)

Solve the coupled major element diffusion equations or the linear diffusion equation for a given domain using finite differences for the discretisation in space and return a solution type variable.

The default time discretisation is based on the ROCK2 method, a stabilized explicit method (Adbdulle and Medovikov, 2001) using OrdinaryDiffEq.

The solution type variable is following the format of OrdinaryDiffEq.jl (see https://docs.sciml.ai/DiffEqDocs/stable/basics/solution/), and can be used to plot the solution, and to extract the solution at a given time. The time of the solution is non-dimensional but can be converted back using the characteristic time (`t_charact` contained in the `Domain` structure).

`path_save` is an optional argument, which can be used to define the path of the HDF5 output file. Default is to nothing.

`solver` is an optional argument, which can be used to define the solver to use for the time discretisation. Default is the ROCK2 method. All other ODE solvers accepted as the ones from DifferentialEquations (see https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).

`p_extra` is an optional argument that can be used to pass additional parameters to the ODEProblem. It is not used by default, but can be useful to write custom callbacks.

All other accepted arguments such as `callback` or `progress` are the same as those of the `solve` function from DifferentialEquations (see https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).
"""
function simulate end

include("simulate_major.jl")
include("simulate_trace.jl")
