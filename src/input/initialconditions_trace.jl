

@kwdef struct IC1DTrace{T1, T2, T3, T4, T5} <: InitialConditions
    C::T1
    Lx::T2
    nx::T3
    Δx::T2
    x::T4
    tfinal::T2
    element::T5
    function IC1DTrace(C::T1, Lx::T2, tfinal::T2, element::String) where {T1 <: AbstractArray{<:Real, 1}, T2 <: Float64}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        # check that element is one of the Rare Earth elements
        elseif !(element in ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"])
            error("Element must be one of the REE: La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu.")
        else
            nx = size(C, 1)
            Δx = Lx / (nx-1)
            x = range(0, length=nx, stop= Lx)

            T3 = typeof(nx)
            T4 = typeof(x)
            T5 = typeof(element)
            new{T1, T2, T3, T4, T5}(C, Lx, nx, Δx, x, tfinal, element)
        end
    end
end


"""
    InitialConditions1DTrace(C::Array{<:Real, 1}, Lx::Unitful.Length, tfinal::Unitful.Time, element::String)

Return a structure containing the initial conditions for a 1D diffusion problem. C needs to be in µg/g. Convert the `Lx` and `tfinal` to µm and Myr respectively.
"""
function InitialConditions1DTrace(;C::AbstractArray{<:Real, 1}, Lx::Unitful.Length, tfinal::Unitful.Time, element::String)
    IC1DTrace(C, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"Myr",tfinal)), element)
end



@kwdef struct Domain1D_trace{T1, T2, T3, T4} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D0::Vector{T2}
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    u0::Vector{T2}
    tfinal_ad::T2
    time_update_ad::T1
    bc_neumann::T3
    function Domain1D_major(IC::IC1DMajor, T::T1, P::T1, time_update::T1, fugacity_O2::T1, bc_neumann::T3) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}, T3 <: Tuple}

        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, Δx, tfinal, Lx, C, element = IC

        T2 = eltype(C)

        u0 = copy(C)

        L_charact = copy(Lx)  # characteristic length
        D_charact = mean(D0)  # characteristic
        t_charact = L_charact^2 / D_charact  # characteristic time

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T4 = typeof(IC)

        new{T1, T2, T3, T4}(IC, T, P, fugacity_O2, time_update, D0, L_charact, D_charact, t_charact, Δxad_, u0, tfinal_ad, time_update_ad, bc_neumann)
    end
end


