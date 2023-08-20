using Unitful
using Unitful: 𝐌, 𝐋, 𝐓, 𝚯
using Parameters

@with_kw struct InitialConditions1D{T <: Real, S <: Real}
    CMg0::Array{Float64, 1}
    CFe0::Array{Float64, 1}
    CMn0::Array{Float64, 1}
    nx::Int
    Lx::T
    Δx::Float64
    x::StepRangeLen
    tfinal::S
    function InitialConditions1D(CMg0, CFe0, CMn0, Lx::T, tfinal::S) where {T <: Real, S <: Real}
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
            new{T,S}(CMg0, CFe0, CMn0, nx, Lx, Δx, x, tfinal)
        end
    end
end

function InitialConditions1D(;CMg0::Array{Float64, 1}, CFe0::Array{Float64, 1}, CMn0::Array{Float64, 1}, Lx::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions1D(CMg0, CFe0, CMn0, ustrip(u"m", Lx), ustrip(u"s", tfinal))
end

@with_kw struct InitialConditions2D{T <: Real, S <: Real}
    CMg0::Array{Float64, 2}
    CFe0::Array{Float64, 2}
    CMn0::Array{Float64, 2}
    nx::Int
    nz::Int
    Lx::T
    Ly::T
    Δx::Float64
    Δy::Float64
    x::StepRangeLen
    y::StepRangeLen
    grid::NamedTuple{(:x, :y), Tuple{Matrix{Float64}, Matrix{Float64}}}
    tfinal::S
    function InitialConditions2D(CMg0, CFe0, CMn0, Lx::T, Ly::T, tfinal::S) where {T <: Real, S <: Real}
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
            new{T,S}(CMg0, CFe0, CMn0, nx, ny, Lx, Ly, Δx, Δy, x, y, grid, tfinal)
        end
    end
end

function InitialConditions2D(;CMg0::Array{Float64, 2}, CFe0::Array{Float64, 2}, CMn0::Array{Float64, 2}, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions2D(CMg0, CFe0, CMn0, ustrip(u"m", Lx), ustrip(u"m", Ly), ustrip(u"s", tfinal))
end

@with_kw struct InitialConditions3D{T <: Real, S <: Real}
    CMg0::Array{Float64, 3}
    CFe0::Array{Float64, 3}
    CMn0::Array{Float64, 3}
    nx::Int
    ny::Int
    nz::Int
    Lx::T
    Ly::T
    Lz::T
    Δx::Float64
    Δy::Float64
    Δz::Float64
    x::StepRangeLen
    y::StepRangeLen
    z::StepRangeLen
    grid::NamedTuple{(:x, :y, :z), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}
    tfinal::S
    function InitialConditions3D(CMg0, CFe0, CMn0, Lx::T, Ly::T, Lz::T, tfinal::S) where {T <: Real, S <: Real}
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
            new{T,S}(CMg0, CFe0, CMn0, nx, ny, nz, Lx, Ly, Lz, Δx, Δy, Δz, x, y, z, grid, tfinal)
        end
    end
end

function InitialConditions3D(;CMg0::Array{Float64, 3}, CFe0::Array{Float64, 3}, CMn0::Array{Float64, 3}, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time)
    InitialConditions3D(CMg0, CFe0, CMn0, ustrip(u"m", Lx), ustrip(u"m", Ly), ustrip(u"m", Lz), ustrip(u"s", tfinal))
end



CMg = ones(100)
CFe = ones(100)
CMn = ones(100)
Lx = 10.0u"m"
Ly = 10.0u"m"
tfinal = 10.0u"s"

test2 = InitialConditions1D(CMg0=CMg, CFe0=CFe, CMn0=CMn,Lx=Lx, tfinal=tfinal)

CMg = ones(100, 100, 100)
CFe = ones(100, 100, 100)
CMn = ones(100, 100, 100)
Lx = 10.0u"m"
Ly = 10.0u"m"
Lz = 10.0u"m"
tfinal = 10.0u"s"

test2 = InitialConditions3D(CMg0=CMg, CFe0=CFe, CMn0=CMn,Lx=Lx,Ly=Ly, Lz=Lz, tfinal=tfinal)

# @with_kw struct Domain{Type <: AbstractRange}

# end

