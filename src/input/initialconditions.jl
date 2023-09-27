using Unitful
using Unitful: 𝐌, 𝐋, 𝐓, 𝚯
using Parameters
using Statistics


@with_kw struct InitialConditions1D{T1, T2, T3, T4}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lx::T2
    nx::T3
    Δx::T2
    x::T4
    tfinal::T2
    function InitialConditions1D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, tfinal::T2) where {T1 <: AbstractArray{<:Real, 1}, T2 <: Float64}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            Δx = Lx / (nx-1)
            x = range(0, length=nx, stop= Lx)

            T3 = typeof(nx)
            T4 = typeof(x)
            new{T1, T2, T3, T4}(CMg0, CFe0, CMn0, Lx, nx, Δx, x, tfinal)
        end
    end
end

@with_kw struct InitialConditionsSpherical{T1, T2, T3, T4}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lr::T2
    nr::T3
    Δr::T2
    r::T4
    tfinal::T2
    function InitialConditionsSpherical(CMg0::T1, CFe0::T1, CMn0::T1, Lr::T2, tfinal::T2) where {T1 <: AbstractArray{<:Real, 1}, T2 <: Float64}
        if Lr <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1)
            error("Initial conditions should have the same size.")
        else
            nr = size(CMg0, 1)
            Δr = Lr / (nr-1)
            radius = range(0, length=nr, stop= Lr)
            T3 = typeof(nr)
            T4 = typeof(radius)
            new{T1, T2, T3, T4}(CMg0, CFe0, CMn0, Lr, nr, Δr, radius, tfinal)
        end
    end
end

@with_kw struct InitialConditions2D{T1, T2, T3, T4}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lx::T2
    Ly::T2
    nx::T3
    ny::T3
    Δx::T2
    Δy::T2
    x::T4
    y::T4
    grt_position::T1
    grt_boundary::T1
    # grid::NamedTuple{(:x, :y), Tuple{AbstractArray{T2, 2}, AbstractArray{T2, 2}}}
    tfinal::T2
    function InitialConditions2D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, tfinal::T2, grt_boundary::T1) where {T1 <: AbstractArray{<:Real, 2}, T2 <: Float64}
        if Lx <= 0 || Ly <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1) || size(CMg0, 2) != size(CFe0, 2) || size(CMg0, 2) != size(CMn0, 2)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            ny = size(CMg0, 2)
            Δx = Lx / (nx-1)
            Δy = Ly / (ny-1)
            x = range(0, length=nx, stop= Lx)
            y = range(0, length=ny, stop= Ly)
            # x first, y second
            # grid = (x=x' .* @ones(ny), y= (@ones(nx))' .* y)

            # define when grt is present
            grt_position = similar(CMg0) .* 0
            # create a binary matrix with 1 where grt is present. 0 otherwise. This depends on if CMg0, CFe0 and CMn0 are equal to 0 or not
            grt_position .= (CMg0 .≠ 0) .& (CFe0 .≠ 0) .& (CMn0 .≠ 0)

            T3 = typeof(nx)
            T4 = typeof(x)

            new{T1, T2, T3, T4}(CMg0, CFe0, CMn0, Lx, Ly, nx, ny, Δx, Δy, x, y, grt_position, grt_boundary, tfinal)
        end
    end
end


@with_kw struct InitialConditions3D{T1 <: AbstractArray{<:Real, 3}, T2 <: Float64}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    nx::Int
    ny::Int
    nz::Int
    Lx::T2
    Ly::T2
    Lz::T2
    Δx::T2
    Δy::T2
    Δz::T2
    x::StepRangeLen
    y::StepRangeLen
    z::StepRangeLen
    tfinal::T2
    function InitialConditions3D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, Lz::T2, tfinal::T2) where {T1 <: AbstractArray{<:Real, 3}, T2 <: Float64}
        if Lx <= 0 || Ly <= 0 || Lz <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1) || size(CMg0, 2) != size(CFe0, 2) || size(CMg0, 2) != size(CMn0, 2) || size(CMg0, 3) != size(CFe0, 3) || size(CMg0, 3) != size(CMn0, 3)
            error("Initial conditions should have the same size.")
        elseif maximum(sum.(CMg0.+ CFe0.+ CMn0)) > 1
            error("Initial conditions should be in mass fraction.")
        else
            nx = size(CMg0, 1)
            ny = size(CMg0, 2)
            nz = size(CMg0, 3)
            Δx = Lx / (nx-1)
            Δy = Ly / (ny-1)
            Δz = Lz / (nz-1)
            x = range(0, length=nx, stop= Lx)
            y = range(0, length=ny, stop= Ly)
            z = range(0, length=nz, stop= Lz)
            new{T1, T2}(CMg0, CFe0, CMn0, nx, ny, nz, Lx, Ly, Lz, Δx, Δy, Δz, x, y, z, tfinal)
        end
    end
end

"""
    InitialConditions1D(CMg0::Array{<:Real, 1}, CFe0::Array{<:Real, 1}, CMn0::Array{<:Real, 1}, Lx::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 1D diffusion problem. CMg0, CFe0 and CMn0 need to be in mass fraction. Convert the `Lx`` and `tfinal`` to µm and Myr respectively.
"""
function InitialConditions1D(CMg0::AbstractArray{<:Real, 1}, CFe0::AbstractArray{<:Real, 1}, CMn0::AbstractArray{<:Real, 1}, Lx::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions1D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"Myr",tfinal)))
end

"""
    InitialConditionsSpherical(CMg0::Array{<:Real, 1}, CFe0::Array{<:Real, 1}, CMn0::Array{<:Real, 1}, Lr::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a spherical diffusion problem. CMg0, CFe0 and CMn0 need to be in mass fraction. Convert `Lr` and `tfinal` to µm and Myr respectively.
"""
function InitialConditionsSpherical(CMg0::AbstractArray{<:Real, 1}, CFe0::AbstractArray{<:Real, 1}, CMn0::AbstractArray{<:Real, 1}, Lr::Unitful.Length, tfinal::Unitful.Time)
    InitialConditionsSpherical(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lr)), convert(Float64,ustrip(u"Myr",tfinal)))
end

"""
    InitialConditions2D(CMg0::Array{<:Real, 2}, CFe0::Array{<:Real, 2}, CMn0::Array{<:Real, 2}, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 2D diffusion problem. CMg0, CFe0 and CMn0 need to be in mass fraction. Convert `Lx`, `Ly` and `tfinal` to µm, µm and Myr respectively.
"""
function InitialConditions2D(CMg0::AbstractArray{<:Real, 2}, CFe0::AbstractArray{<:Real, 2}, CMn0::AbstractArray{<:Real, 2}, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time; grt_boundary::AbstractArray{<:Real, 2}=similar(CMg0) .* 0.0)
    InitialConditions2D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"µm", Ly)), convert(Float64,ustrip(u"Myr", tfinal)), grt_boundary)
end

"""
    InitialConditions3D(CMg0::Array{<:Real, 3}, CFe0::Array{<:Real, 3}, CMn0::Array{<:Real, 3}, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 3D diffusion problem. CMg0, CFe0 and CMn0 need to be in mass fraction. Convert `Lx`, `Ly`, `Lz` and `tfinal` to µm, µm, µm and Myr respectively.
"""
function InitialConditions3D(CMg0::AbstractArray{<:Real, 3}, CFe0::AbstractArray{<:Real, 3}, CMn0::AbstractArray{<:Real, 3}, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions3D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"µm", Ly)), convert(Float64,ustrip(u"µm", Lz)), convert(Float64,ustrip(u"Myr", tfinal)))
end


function D_ini!(D0,T,P)
    R = 8.314462618  # constant of gas in J mol−1 K−1

    # Magnesium
    D0Mg = 1.1 * 1e-3 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Mg = 67997 * 4.1855 # activation energy at 1 bar in J / mol
    ΔV⁺Mg = 5.3  # activation volume in cm3 / mol

    # Iron
    D0Fe = 6.4 * 1e-4 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Fe = 65824 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Fe = 5.6  # activation volume in cm3 / mol

    # Manganese
    D0Mn = 5.1 * 1e-4 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Mn = 60569 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Mn = 6.0  # activation volume in cm3 / mol

    DMg = D0Mg * exp(- (Eₐ_Mg + (100 * (P-0.001) * ΔV⁺Mg)) / (R * (T+273.15)))  # in µm2 / s
    DFe = D0Fe * exp(- (Eₐ_Fe + (100 * (P-0.001) * ΔV⁺Fe)) / (R * (T+273.15)))  # in µm2 / s
    DMn = D0Mn * exp(- (Eₐ_Mn + (100 * (P-0.001) * ΔV⁺Mn)) / (R * (T+273.15)))  # in µm2 / s
    DCa = 0.5 * DFe


    D0 .= [DMg, DFe, DMn, DCa] .* (365.25 * 24 * 3600 * 1e6)  # in years
end

@with_kw struct Domain1D{T1 <: Union{AbstractArray{Float64, 1}, Float64}, T2 <: Float64, T3 <: Tuple}
    IC::InitialConditions1D
    T::T1
    P::T1
    time_update::T1
    D0::Vector{T2}
    D::NamedTuple{(:DMgMg, :DMgFe, :DMgMn, :DFeMg, :DFeFe, :DFeMn, :DMnMg, :DMnFe, :DMnMn),
                  NTuple{9, Vector{T2}}}  # tensor of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    u0::Matrix{T2}
    tfinal_ad::T2
    time_update_ad::T1
    bc_neumann::T3
    function Domain1D(IC::InitialConditions1D, T::T1, P::T1, time_update::T1, bc_neumann::T2) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}, T2 <: Tuple}

        #check that T, P and time_update have the same size
        if size(T, 1) != size(P, 1) || size(T, 1) != size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, Δx, tfinal, Lx, CMg0, CFe0, CMn0 = IC

        D0 = zeros(Float64, 4)
        D_ini!(D0, T[1], P[1])  # compute initial diffusion coefficients

        D = (DMgMg = zeros(nx), DMgFe = zeros(nx), DMgMn = zeros(nx), DFeMg = zeros(nx), DFeFe = zeros(nx), DFeMn = zeros(nx), DMnMg = zeros(nx), DMnFe = zeros(nx), DMnMn = zeros(nx))  # tensor of interdiffusion coefficients

        u0 = similar(CMg0, (nx, 3))
        u0[:,1] .= CMg0
        u0[:,2] .= CFe0
        u0[:,3] .= CMn0

        L_charact = copy(Lx)  # characteristic length
        D_charact = mean(D0)  # characteristic
        t_charact = L_charact^2 / D_charact  # characteristic time

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update
        new{T1, Float64, T2}(IC, T, P, time_update, D0, D, L_charact, D_charact, t_charact, Δxad_, u0, tfinal_ad, time_update_ad, bc_neumann)
    end
end

@with_kw struct DomainSpherical{T1 <: Union{AbstractArray{Float64, 1}, Float64}, T2 <: Float64}
    IC::InitialConditionsSpherical
    T::T1
    P::T1
    time_update::T1
    D0::Vector{T2}
    D::NamedTuple{(:DMgMg, :DMgFe, :DMgMn, :DFeMg, :DFeFe, :DFeMn, :DMnMg, :DMnFe, :DMnMn),
                  NTuple{9, Vector{T2}}}  # tensor of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δrad::T2
    Δrad_::T2
    r_ad::Vector{Float64}
    u0::Matrix{T2}
    tfinal_ad::T2
    time_update_ad::T1
    function DomainSpherical(IC::InitialConditionsSpherical, T::T1, P::T1, time_update::T1) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}

        #check that T, P and time_update have the same size
        if size(T, 1) != size(P, 1) || size(T, 1) != size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nr, Δr, r, tfinal, Lr, CMg0, CFe0, CMn0 = IC

        D0::Vector{Float64} = zeros(Float64, 4)
        D_ini!(D0, T[1], P[1])  # compute initial diffusion coefficients

        D = (DMgMg = zeros(nr), DMgFe = zeros(nr), DMgMn = zeros(nr), DFeMg = zeros(nr), DFeFe = zeros(nr), DFeMn = zeros(nr), DMnMg = zeros(nr), DMnFe = zeros(nr), DMnMn = zeros(nr))  # tensor of interdiffusion coefficients

        u0::Matrix{typeof(CMg0[1])} = similar(CMg0, (nr, 3))
        u0[:,1] .= CMg0
        u0[:,2] .= CFe0
        u0[:,3] .= CMn0

        L_charact = copy(Lr)  # characteristic length
        D_charact = mean(D0)  # characteristic
        t_charact = L_charact^2 / D_charact  # characteristic time

        Δrad = Δr / L_charact  # nondimensionalised Δr
        Δrad_ = 1 / Δrad  # inverse of nondimensionalised Δr
        r_ad::Vector{Float64} = r ./ L_charact  # nondimensionalised radius
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update
        new{T1, Float64}(IC, T, P, time_update, D0, D, L_charact, D_charact, t_charact, Δrad, Δrad_, r_ad, u0, tfinal_ad, time_update_ad)
    end
end

@with_kw struct Domain2D{T1, T2, T3, T4, T5, T6}
    IC::InitialConditions2D
    T::T1
    P::T1
    time_update::T1
    D0::T4
    D::T5  # tensor of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    Δyad_::T2
    u0::T6
    tfinal_ad::T2
    function Domain2D(IC::InitialConditions2D, T::T1, P::T1, time_update::T1) where {T1 <: Union{Float64, AbstractArray{Float64, 2}}}
        @unpack nx, ny, Δx, Δy, tfinal, Lx, CMg0, CFe0, CMn0 = IC
        similar(CMg0, (nx,ny))
        D0 = similar(CMg0, (4)) .* 0
        D_ini!(D0, T, P)  # compute initial diffusion coefficients

        D = (DMgMg = similar(CMg0, (nx,ny)) .* 0, DMgFe = similar(CMg0, (nx,ny)) .* 0, DMgMn = similar(CMg0, (nx,ny)) .* 0, DFeMg = similar(CMg0, (nx,ny)) .* 0, DFeFe = similar(CMg0, (nx,ny)) .* 0, DFeMn = similar(CMg0, (nx,ny)) .* 0, DMnMg = similar(CMg0, (nx,ny)) .* 0, DMnFe = similar(CMg0, (nx,ny)) .* 0, DMnMn = similar(CMg0, (nx,ny)) .* 0)  # tensor of interdiffusion coefficients

        u0 = similar(CMg0, (nx, ny, 3)) .* 0
        u0[:, :, 1] .= CMg0
        u0[:, :, 2] .= CFe0
        u0[:, :, 3] .= CMn0

        L_charact = copy(Lx)  # characteristic length
        D_charact = mean(D0)  # characteristic
        t_charact = L_charact^2 / D_charact  # characteristic time
        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        Δyad_ = 1 / (Δy / L_charact)  # inverse of nondimensionalised Δy
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update = time_update / t_charact  # nondimensionalised time update

        T2 = typeof(t_charact)
        T3 = typeof(CMg0)
        T4 = typeof(D0)
        T5 = typeof(D)
        T6 = typeof(u0)

        new{T1, T2, T3, T4, T5, T6}(IC, T, P, time_update, D0, D, L_charact, D_charact, t_charact, Δxad_, Δyad_, u0, tfinal_ad)
    end
end

@with_kw struct Domain3D{T1 <: Union{Array{Float64, 3}, Float64}, T2 <: Float64}
    IC::InitialConditions3D
    T::T1
    P::T1
    time_update::T1
    D0::Vector{T2}
    D::NamedTuple{(:DMgMg, :DMgFe, :DMgMn, :DFeMg, :DFeFe, :DFeMn, :DMnMg, :DMnFe, :DMnMn),
                  Tuple{Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}, Array{T2, 3}}}  # tensor of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    Δyad_::T2
    Δzad_::T2
    u0::Array{T2, 4}
    tfinal_ad::T2
    function Domain3D(IC::InitialConditions3D, T::T1, P::T1, time_update::T1) where {T1 <: Union{Float64, Array{Float64, 3}}}
        @unpack nx, ny, nz, Δx, Δy, Δz, tfinal, Lx, CMg0, CFe0, CMn0 = IC

        D0 = @zeros(4)
        D_ini!(D0, T, P)  # compute initial diffusion coefficients

        D = (DMgMg = zeros(nx, ny, nz), DMgFe = zeros(nx, ny, nz), DMgMn = zeros(nx, ny, nz), DFeMg = zeros(nx, ny, nz), DFeFe = zeros(nx, ny, nz), DFeMn = zeros(nx, ny, nz), DMnMg = zeros(nx, ny, nz), DMnFe = zeros(nx, ny, nz), DMnMn = zeros(nx, ny, nz))  # tensor of interdiffusion coefficients


        u0 = similar(CMg0, (nx, ny, nz, 3))
        u0[:, :, :, 1] .= CMg0
        u0[:, :, :, 2] .= CFe0
        u0[:, :, :, 3] .= CMn0

        L_charact = copy(Lx)  # characteristic length
        D_charact = mean(D0)  # characteristic
        t_charact = L_charact^2 / D_charact  # characteristic time

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        Δyad_ = 1 / (Δy / L_charact)  # inverse of nondimensionalised Δy
        Δzad_ = 1 / (Δz / L_charact)  # inverse of nondimensionalised Δz
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update = time_update / t_charact  # nondimensionalised time update
        new{T1, Float64}(IC, T, P, time_update, D0, D, L_charact, D_charact, t_charact, Δxad_, Δyad_, Δzad_, u0, tfinal_ad)
    end
end

"""
    Domain(IC, T, P, time_update=0u"Myr")

Return a struct containing the struct IC, the PT conditions `T` and `P`  and the nondimensionalised parameters of the problem. `time_update` can be used to update the PT conditions during the simulation. If so, `T` and `P` have to be of the same size as `time_update` and correspond to the values of PT to update to.

"""
function Domain end

"""
    Domain(IC::InitialConditions1D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr"; bc_neumann::Tuple=(false, false))

When applied to 1D initial conditions, define corresponding 1D domain. `bc_neumann` can be used to define Neumann boundary conditions on the left or right of the domain if set to true.
"""
function Domain(IC::InitialConditions1D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr"; bc_neumann::Tuple=(false, false))
    Domain1D(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), bc_neumann)
end

"""
   Domain(IC::InitialConditionsSpherical, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to spherical initial conditions, define corresponding spherical domain.
"""
function Domain(IC::InitialConditionsSpherical, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")
    DomainSpherical(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)))
end

"""
    Domain(IC::InitialConditions2D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 2D initial conditions, define corresponding 2D domain.
"""
function Domain(IC::InitialConditions2D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")
    Domain2D(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)))
end

"""
    Domain(IC::InitialConditions3D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 3D initial conditions, define corresponding 3D domain.
"""
function Domain(IC::InitialConditions3D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")
    Domain3D(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)))
end