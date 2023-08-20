using Unitful
using Unitful: 𝐌, 𝐋, 𝐓, 𝚯
using Parameters
using Statistics


@with_kw struct InitialConditions1D{T1<:Array{Float64, 1} , T2 <: Float64}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    nx::Int
    Lx::T2
    Δx::T2
    x::StepRangeLen
    tfinal::T2
    function InitialConditions1D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, tfinal::T2) where {T1 <: Array{Float64, 1}, T2 <: Float64}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            Δx = Lx / (nx)
            x = range(Δx/2, length=nx, stop= Lx-Δx/2)
            new{T1, T2}(CMg0, CFe0, CMn0, nx, Lx, Δx, x, tfinal)
        end
    end
end

function InitialConditions1D(;CMg0::Array{Float64, 1}, CFe0::Array{Float64, 1}, CMn0::Array{Float64, 1}, Lx::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions1D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"m", Lx)), convert(Float64,ustrip(u"Myr",tfinal)))
end

@with_kw struct InitialConditions2D{T1 <: Array{Float64, 2}, T2 <: Float64}
    CMg0::T1
    CFe0::T1
    CMn0::T1
    nx::Int
    ny::Int
    Lx::T2
    Ly::T2
    Δx::T2
    Δy::T2
    x::StepRangeLen
    y::StepRangeLen
    grid::NamedTuple{(:x, :y), Tuple{Matrix{Float64}, Matrix{Float64}}}
    tfinal::T2
    function InitialConditions2D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, tfinal::T2) where {T1 <: Array{Float64, 2}, T2 <: Float64}
        if Lx <= 0 || Ly <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1) || size(CMg0, 2) != size(CFe0, 2) || size(CMg0, 2) != size(CMn0, 2)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            ny = size(CMg0, 2)
            Δx = Lx / (nx)
            Δy = Ly / (ny)
            x = range(Δx/2, length=nx, stop= Lx-Δx/2)
            y = range(Δy/2, length=ny, stop= Ly-Δy/2)
            # x first, y second
            grid = (x=x' .* ones(ny), y=ones(nx)' .* y)
            new{T1, T2}(CMg0, CFe0, CMn0, nx, ny, Lx, Ly, Δx, Δy, x, y, grid, tfinal)
        end
    end
end

function InitialConditions2D(;CMg0::Array{Float64, 2}, CFe0::Array{Float64, 2}, CMn0::Array{Float64, 2}, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions2D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"m", Lx)), convert(Float64,ustrip(u"m", Ly)), convert(Float64,ustrip(u"Myr", tfinal)))
end


@with_kw struct InitialConditions3D{T1 <: Array{Float64, 3}, T2 <: Float64}
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
    grid::NamedTuple{(:x, :y, :z), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}
    tfinal::T2
    function InitialConditions3D(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, Lz::T2, tfinal::T2) where {T1 <: Array{Float64, 3}, T2 <: Float64}
        if Lx <= 0 || Ly <= 0 || Lz <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) != size(CFe0, 1) || size(CMg0, 1) != size(CMn0, 1) || size(CMg0, 2) != size(CFe0, 2) || size(CMg0, 2) != size(CMn0, 2) || size(CMg0, 3) != size(CFe0, 3) || size(CMg0, 3) != size(CMn0, 3)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            ny = size(CMg0, 2)
            nz = size(CMg0, 3)
            Δx = Lx / (nx)
            Δy = Ly / (ny)
            Δz = Lz / (nz)
            x = range(Δx/2, length=nx, stop= Lx-Δx/2)
            y = range(Δy/2, length=ny, stop= Ly-Δy/2)
            z = range(Δz/2, length=nz, stop= Lz-Δz/2)
            # x first, y second, z third
            grid = (x=x' .* ones(ny, nz), y=ones(nx)' .* y .* ones(nz), z=ones(nx, ny) .* z)
            new{T1, T2}(CMg0, CFe0, CMn0, nx, ny, nz, Lx, Ly, Lz, Δx, Δy, Δz, x, y, z, grid, tfinal)
        end
    end
end

function InitialConditions3D(;CMg0::Array{Float64, 3}, CFe0::Array{Float64, 3}, CMn0::Array{Float64, 3}, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions3D(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"m", Lx)), convert(Float64,ustrip(u"m", Ly)), convert(Float64,ustrip(u"m", Lz)), convert(Float64,ustrip(u"Myr", tfinal)))
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

@with_kw struct Domain1D{T1 <: Float64}
    IC::InitialConditions1D
    T::T1
    P::T1
    time_update::T1
    D0::Vector{T1}
    L_charact::T1
    D_charact::T1
    t_charact::T1
    Δxad_::T1
    tfinal_ad::T1
    function Domain1D(IC::InitialConditions1D, T::T1, P::T1, time_update::T1) where {T1 <: Float64}
        D0 = zeros(Float64, 4)
        @show T
        D_ini!(D0, T, P)
        @show D0
        L_charact = copy(IC.Lx)
        D_charact = mean(D0)
        t_charact = L_charact^2 / D_charact
        Δxad_ = 1 / (IC.Δx / L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = IC.tfinal / t_charact  # nondimensionalised total time
        time_update = time_update / t_charact  # nondimensionalised time update
        new{T1}(IC, T, P, time_update, D0, L_charact, D_charact, t_charact, Δxad_, tfinal_ad)
    end
end

function Domain1D(;IC::InitialConditions1D, T::Unitful.Temperature, P::Unitful.Pressure, time_update::Unitful.Time=0u"Myr")
    Domain1D(IC, convert(Float64,ustrip(u"°C", T)), convert(Float64,ustrip(u"kbar", P)), convert(Float64,ustrip(u"Myr", time_update)))
end
