# [Updating pressure and temperature conditions](@id PT_callback)

When modelling diffusion in garnet, it may be relevant to update the pressure and temperature (PT) conditions during the simulation.

To do this, we can use a callback function. A callback function is a function that can be called in our solver when a certain condition is met, i.e when a certain time is reached. It is based on the library [DiffEqCallbacks](https://docs.sciml.ai/DiffEqCallbacks/stable/) from the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) ecosystem.