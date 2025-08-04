using Unitful
using Unitful:𝐋,𝐓
using Statistics

# define diffusion coefficient as a derived dimension
@derived_dimension 𝓓 𝐋^2/𝐓

abstract type InitialConditions end
abstract type Domain end

"""
    Domain(IC, T, P, time_update=0u"Myr")

Return a struct containing the struct `IC`, the PT conditions `T` and `P`  and the nondimensionalised parameters of the problem. `time_update` can be used to update the PT conditions during the simulation. If so, `T` and `P` have to be of the same size as `time_update` and correspond to the values of PT to update to.

"""
function Domain end

include("initialconditions_major.jl")
# include("initialconditions_trace.jl")
