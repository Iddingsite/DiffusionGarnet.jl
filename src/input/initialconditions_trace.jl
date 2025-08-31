
@kwdef struct InitialConditions1DTrace{T1, T2, T_float, T_Int, T_Vec} <: InitialConditions
    C::T1
    D::T2
    Lx::T_float
    nx::T_Int
    Δx::T_Vec
    x::T_Vec
    tfinal::T_float
    function InitialConditions1DTrace(C::T1, D::T2, Lx::T3, x::ArrayX, tfinal::T3) where {T1 <: AbstractArray{<:Real, 1}, T2 <: AbstractChemicalDiffusion, T3 <: Float64, ArrayX <: AbstractArray{<:Real, 1}}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            nx = size(C, 1)
            Δx = diff(x)

            T4 = typeof(nx)

            new{T1, T2, T3, T4, ArrayX}(C, D, Lx, nx, Δx, x, tfinal)
        end
    end
end


"""
    IC1DTrace(;C::Array{<:Real, 1}, D, Lx::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 1D diffusion problem. C needs to be in µg/g. Convert the `Lx` and `tfinal` to cm and Myr respectively.
"""
function IC1DTrace(;C::AbstractArray{<:Real, 1},
                    D,
                    Lx::Unitful.Length,
                    x::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 1), stop=Lx),
                    tfinal::Unitful.Time
                    )

    InitialConditions1DTrace(C, D, convert(Float64,ustrip(u"µm", Lx)), ustrip.(u"µm", x), convert(Float64,ustrip(u"Myr",tfinal)))
end


@kwdef struct Domain1DTrace{T1, T2, T3, T4, T_Vec} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D::T2
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T_Vec
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
        D =  ustrip(u"µm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar))  # diffusion coefficient in µm^2/Myr

        L_charact = copy(Lx)  # characteristic length
        t_charact = 1.0  # characteristic time
        D_charact = L_charact^2 / t_charact  # characteristic diffusion coefficient

        Δxad_ = 1 ./ (Δx ./ L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T4 = typeof(IC)

        T_vec =  typeof(Δxad_)

        new{T1, T2, T3, T4, T_vec}(IC, T, P, fugacity_O2, time_update, D, L_charact, D_charact, t_charact, Δxad_, u0, tfinal_ad, time_update_ad, bc_neumann)
    end
end


"""
    Domain(IC::InitialConditions1DTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 1D initial conditions, define corresponding 1D domain for trace element diffusion. `bc_neumann` can be used to define Neumann boundary conditions on the left or right side of the domain if set to true.
"""
function Domain(IC::InitialConditions1DTrace,
                T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}},
                P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}},
                time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr",
                fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa";
                bc_neumann::Tuple=(false, false)
                )

    Domain1DTrace(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), bc_neumann)
end

@kwdef struct InitialConditionsSphericalTrace{T1, T2, T_float, T_Int, T_Vec} <: InitialConditions
    C::T1
    D::T2
    Lr::T_float
    nr::T_Int
    Δr::T_Vec
    r::T_Vec
    tfinal::T_float
    function InitialConditionsSphericalTrace(C::T1, D::T2, Lr::T3, r::ArrayR, tfinal::T3) where {T1 <: AbstractArray{<:Real, 1}, T2 <: AbstractChemicalDiffusion, T3 <: Float64, ArrayR <: AbstractArray{<:Real, 1}}
        if Lr <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            nr = size(C, 1)
            Δr = diff(r)
            T4 = typeof(nr)
            new{T1, T2, T3, T4, ArrayR}(C, D, Lr, nr, Δr, r, tfinal)
        end
    end
end

"""
    ICSphTrace(;C::Array{<:Real, 1}, D, Lr::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a spherical diffusion problem. C needs to be in µg/g. Convert `Lr` and `tfinal` to µm and Myr respectively.
"""
function ICSphTrace(;C::AbstractArray{<:Real, 1},
                    D,
                    Lr::Unitful.Length,
                    r::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 1), stop=Lr),
                    tfinal::Unitful.Time
                    )
    InitialConditionsSphericalTrace(C, D, convert(Float64,ustrip(u"µm", Lr)), ustrip.(u"µm", r), convert(Float64,ustrip(u"Myr",tfinal)))
end

@kwdef struct DomainSphericalTrace{T1, T2, T3, T_Vec} <: Domain
    IC::T3
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D::T2
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δr_ad::T_Vec
    Δr_ad_::T_Vec
    r_ad::T_Vec
    u0::Vector{T2}
    tfinal_ad::T2
    time_update_ad::T1
    function DomainSphericalTrace(IC::InitialConditionsSphericalTrace, T::T1, P::T1, time_update::T1, fugacity_O2::T1) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}
        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nr, Δr, r, tfinal, Lr, C = IC

        T2 = eltype(C)
        u0 = copy(C)

        T_K = (T+273.15) * 1u"K"
        P_kbar = P * 1u"kbar"
        D = ustrip(u"µm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar))  # diffusion coefficient in µm^2/Myr

        L_charact = copy(Lr)  # characteristic length
        t_charact = 1.0  # characteristic time
        D_charact = L_charact^2 / t_charact  # characteristic diffusion coefficient

        Δr_ad = Δr ./ L_charact  # nondimensionalised Δr
        Δr_ad_ = 1 ./ Δr_ad      # inverse of nondimensionalised Δr
        r_ad = r ./ L_charact    # nondimensionalised radius
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T3 = typeof(IC)
        T_vec = typeof(Δr_ad)

        new{T1, T2, T3, T_vec}(IC, T, P, fugacity_O2, time_update, D, L_charact, D_charact, t_charact, Δr_ad, Δr_ad_, r_ad, u0, tfinal_ad, time_update_ad)
    end
end

"""
    Domain(IC::InitialConditionsSphericalTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to spherical initial conditions, define corresponding spherical domain for trace element diffusion.
"""
function Domain(IC::InitialConditionsSphericalTrace,
                T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}},
                P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}},
                time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr",
                fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa";
                )

    DomainSphericalTrace(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)))
end


@kwdef struct InitialConditions2DTrace{T1, T2, T_float, T_Int, T_VecX, T_VecY} <: InitialConditions
    C::T1
    D::T2
    Lx::T_float
    Ly::T_float
    nx::T_Int
    ny::T_Int
    Δx::T_VecX
    Δy::T_VecY
    x::T_VecX
    y::T_VecY
    tfinal::T_float
    function InitialConditions2DTrace(C::T1, D::T2, Lx::T3, Ly::T3, x::ArrayX, y::ArrayY, tfinal::T3) where {T1 <: AbstractArray{<:Real, 2}, T2 <: AbstractChemicalDiffusion, T3 <: Float64, ArrayX <: AbstractArray{<:Real, 1}, ArrayY <: AbstractArray{<:Real, 1}}
        if Lx <= 0 || Ly <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            nx = size(C, 1)
            ny = size(C, 2)
            Δx = diff(x)
            Δy = diff(y)
            T4 = typeof(nx)
            T5 = typeof(ny)
            new{T1, T2, T3, T4, ArrayX, ArrayY}(C, D, Lx, Ly, nx, ny, Δx, Δy, x, y, tfinal)
        end
    end
end

"""
    IC2DTrace(;C::Array{<:Real, 2}, D, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 2D diffusion problem. C needs to be in µg/g. Convert `Lx`, `Ly` and `tfinal` to µm, µm and Myr respectively.
"""
function IC2DTrace(;C::AbstractArray{<:Real, 2},
                   D,
                   Lx::Unitful.Length,
                   Ly::Unitful.Length,
                   x::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 1), stop=Lx),
                   y::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 2), stop=Ly),
                   tfinal::Unitful.Time
                   )
    InitialConditions2DTrace(C, D, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"µm", Ly)), ustrip.(u"µm", x), ustrip.(u"µm", y), convert(Float64,ustrip(u"Myr",tfinal)))
end

@kwdef struct Domain2DTrace{T1, T2, T3, T4, T_VecX, T_VecY} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D::T2
    L_charact_x::T2
    L_charact_y::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T_VecX
    Δyad_::T_VecY
    u0::Vector{T2}
    tfinal_ad::T2
    time_update_ad::T1
    function Domain2DTrace(IC::InitialConditions2DTrace, T::T1, P::T1, time_update::T1, fugacity_O2::T1) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}
        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, ny, Δx, Δy, x, y, tfinal, Lx, Ly, C = IC

        T2 = eltype(C)
        u0 = copy(C)

        T_K = (T+273.15) * 1u"K"
        P_kbar = P * 1u"kbar"
        D = ustrip(u"µm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar))  # diffusion coefficient in µm^2/Myr

        L_charact_x = copy(Lx)
        L_charact_y = copy(Ly)
        t_charact = 1.0
        D_charact = (L_charact_x^2 + L_charact_y^2) / (2 * t_charact)

        Δxad_ = 1 ./ (Δx ./ L_charact_x)
        Δyad_ = 1 ./ (Δy ./ L_charact_y)
        tfinal_ad = tfinal / t_charact
        time_update_ad = time_update ./ t_charact

        T4 = typeof(IC)
        T_vecx = typeof(Δxad_)
        T_vecy = typeof(Δyad_)

        new{T1, T2, T3, T4, T_vecx, T_vecy}(IC, T, P, fugacity_O2, time_update, D, L_charact_x, L_charact_y, D_charact, t_charact, Δxad_, Δyad_, u0, tfinal_ad, time_update_ad)
    end
end

"""
    Domain(IC::InitialConditions2DTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 2D initial conditions, define corresponding 2D domain for trace element diffusion.
"""
function Domain(IC::InitialConditions2DTrace,
                T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}},
                P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}},
                time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr",
                fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa";
                )

    Domain2DTrace(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)))
end

@kwdef struct InitialConditions3DTrace{T1, T2, T_float, T_Int, T_VecX, T_VecY, T_VecZ} <: InitialConditions
    C::T1
    D::T2
    Lx::T_float
    Ly::T_float
    Lz::T_float
    nx::T_Int
    ny::T_Int
    nz::T_Int
    Δx::T_VecX
    Δy::T_VecY
    Δz::T_VecZ
    x::T_VecX
    y::T_VecY
    z::T_VecZ
    tfinal::T_float
    function InitialConditions3DTrace(C::T1, D::T2, Lx::T3, Ly::T3, Lz::T3, x::ArrayX, y::ArrayY, z::ArrayZ, tfinal::T3) where {T1 <: AbstractArray{<:Real, 3}, T2 <: AbstractChemicalDiffusion, T3 <: Float64, ArrayX <: AbstractArray{<:Real, 1}, ArrayY <: AbstractArray{<:Real, 1}, ArrayZ <: AbstractArray{<:Real, 1}}
        if Lx <= 0 || Ly <= 0 || Lz <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            nx = size(C, 1)
            ny = size(C, 2)
            nz = size(C, 3)
            Δx = diff(x)
            Δy = diff(y)
            Δz = diff(z)
            T4 = typeof(nx)
            T5 = typeof(ny)
            T6 = typeof(nz)
            new{T1, T2, T3, T4, ArrayX, ArrayY, ArrayZ}(C, D, Lx, Ly, Lz, nx, ny, nz, Δx, Δy, Δz, x, y, z, tfinal)
        end
    end
end

"""
    IC3DTrace(;C::Array{<:Real, 3}, D, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 3D diffusion problem. C needs to be in µg/g. Convert `Lx`, `Ly`, `Lz` and `tfinal` to µm, µm, µm and Myr respectively.
"""
function IC3DTrace(;C::AbstractArray{<:Real, 3},
                   D,
                   Lx::Unitful.Length,
                   Ly::Unitful.Length,
                   Lz::Unitful.Length,
                   x::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 1), stop=Lx),
                   y::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 2), stop=Ly),
                   z::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(C, 3), stop=Lz),
                   tfinal::Unitful.Time
                   )
    InitialConditions3DTrace(C, D,
        convert(Float64,ustrip(u"µm", Lx)),
        convert(Float64,ustrip(u"µm", Ly)),
        convert(Float64,ustrip(u"µm", Lz)),
        ustrip.(u"µm", x),
        ustrip.(u"µm", y),
        ustrip.(u"µm", z),
        convert(Float64,ustrip(u"Myr",tfinal)))
end

@kwdef struct Domain3DTrace{T1, T2, T3, T4, T_VecX, T_VecY, T_VecZ} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    D::T2
    L_charact_x::T2
    L_charact_y::T2
    L_charact_z::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T_VecX
    Δyad_::T_VecY
    Δzad_::T_VecZ
    u0::Vector{T2}
    tfinal_ad::T2
    time_update_ad::T1
    function Domain3DTrace(IC::InitialConditions3DTrace, T::T1, P::T1, time_update::T1, fugacity_O2::T1) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}
        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, ny, nz, Δx, Δy, Δz, x, y, z, tfinal, Lx, Ly, Lz, C = IC

        T2 = eltype(C)
        u0 = copy(C)

        T_K = (T+273.15) * 1u"K"
        P_kbar = P * 1u"kbar"
        D = ustrip(u"µm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar))  # diffusion coefficient in µm^2/Myr

        L_charact_x = copy(Lx)
        L_charact_y = copy(Ly)
        L_charact_z = copy(Lz)
        t_charact = 1.0
        D_charact = (L_charact_x^2 + L_charact_y^2 + L_charact_z^2) / (3 * t_charact)

        Δxad_ = 1 ./ (Δx ./ L_charact_x)
        Δyad_ = 1 ./ (Δy ./ L_charact_y)
        Δzad_ = 1 ./ (Δz ./ L_charact_z)
        tfinal_ad = tfinal / t_charact
        time_update_ad = time_update ./ t_charact

        T4 = typeof(IC)
        T_vecx = typeof(Δxad_)
        T_vecy = typeof(Δyad_)
        T_vecz = typeof(Δzad_)

        new{T1, T2, T3, T4, T_vecx, T_vecy, T_vecz}(IC, T, P, fugacity_O2, time_update, D, L_charact_x, L_charact_y, L_charact_z, D_charact, t_charact, Δxad_, Δyad_, Δzad_, u0, tfinal_ad, time_update_ad)
    end
end

"""
    Domain(IC::InitialConditions3DTrace, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 3D initial conditions, define corresponding 3D domain for trace element diffusion.
"""
function Domain(IC::InitialConditions3DTrace,
                T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}},
                P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}},
                time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr",
                fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa";
                )

    Domain3DTrace(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)))
end

