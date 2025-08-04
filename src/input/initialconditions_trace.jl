
@kwdef struct InitialConditions1DTrace{T1, T2, T3, T4, T5} <: InitialConditions
    C::T1
    D::T2
    Lx::T3
    nx::T4
    Δx::T3
    x::T5
    tfinal::T3
    function InitialConditions1DTrace(C::T1, D::T2, Lx::T3, tfinal::T3) where {T1 <: AbstractArray{<:Real, 1}, T2 <: AbstractChemicalDiffusion, T3 <: Float64}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            nx = size(C, 1)
            Δx = Lx / (nx-1)
            x = range(0, length=nx, stop= Lx)

            T4 = typeof(nx)
            T5 = typeof(x)
            new{T1, T2, T3, T4, T5}(C, D, Lx, nx, Δx, x, tfinal)
        end
    end
end


"""
    IC1DTrace(;C::Array{<:Real, 1}, D, Lx::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 1D diffusion problem. C needs to be in µg/g. Convert the `Lx` and `tfinal` to cm and Myr respectively.
"""
function IC1DTrace(;C::AbstractArray{<:Real, 1}, D, Lx::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions1DTrace(C, D, convert(Float64,ustrip(u"cm", Lx)), convert(Float64,ustrip(u"Myr",tfinal)))
end


@kwdef struct Domain1DTrace{T1, T2, T3, T4} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D::T2
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    u0::Vector{T2}
    tfinal_ad::T2
    time_update_ad::T1
    bc_neumann::T3
    function Domain1DTrace(IC::InitialConditions1DTrace, T::T1, P::T1, time_update::T1, fugacity_O2::T1, bc_neumann::T3) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}, T3 <: Tuple}

        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, Δx, tfinal, Lx, C = IC

        T2 = eltype(C)

        u0 = copy(C)

        T_K = (T+273.15) * 1u"K"
        P_kbar = P * 1u"kbar"
        D =  ustrip(u"cm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar))  # diffusion coefficient in µm^2/Myr

        L_charact = copy(Lx)  # characteristic length
        D_charact = D  # characteristic diffusion coefficient
        t_charact = L_charact^2 / D_charact  # characteristic time

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T4 = typeof(IC)

        new{T1, T2, T3, T4}(IC, T, P, fugacity_O2, time_update, D, L_charact, D_charact, t_charact, Δxad_, u0, tfinal_ad, time_update_ad, bc_neumann)
    end
end


"""
    Domain(IC::InitialConditions1DTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 1D initial conditions, define corresponding 1D domain for trace element diffusion. `bc_neumann` can be used to define Neumann boundary conditions on the left or right side of the domain if set to true.
"""
function Domain(IC::InitialConditions1DTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr", fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa"; bc_neumann::Tuple=(false, false))
    Domain1DTrace(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), bc_neumann)
end