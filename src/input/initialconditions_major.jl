
@kwdef struct InitialConditions1DMajor{T1, T2, T3, T4} <: InitialConditions
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lx::T2
    nx::T3
    Δx::T4
    x::T4
    tfinal::T2
    function InitialConditions1DMajor(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, x::ArrayX, tfinal::T2) where {T1 <: AbstractArray{<:Real, 1}, T2 <: Float64, ArrayX <: AbstractArray{<:Real, 1}}
        if Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) ≠ size(CFe0, 1) || size(CMg0, 1) ≠ size(CMn0, 1)
            error("Initial conditions should have the same size.")
        else
            nx = size(CMg0, 1)
            Δx = diff(x)

            T3 = typeof(nx)

            new{T1, T2, T3, ArrayX}(CMg0, CFe0, CMn0, Lx, nx, Δx, x, tfinal)
        end
    end
end

@kwdef struct InitialConditionsSpherical{T1, T2, T3, T4} <: InitialConditions
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lr::T2
    nr::T3
    Δr::T4
    r::T4
    tfinal::T2
    function InitialConditionsSpherical(CMg0::T1, CFe0::T1, CMn0::T1, Lr::T2, r::ArrayR, tfinal::T2) where {T1 <: AbstractArray{<:Real, 1}, T2 <: Float64, ArrayR <: AbstractArray{<:Real, 1}}
        if Lr <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) ≠ size(CFe0, 1) || size(CMg0, 1) ≠ size(CMn0, 1)
            error("Initial conditions should have the same size.")
        else
            nr = size(CMg0, 1)
            Δr = diff(r)

            T3 = typeof(nr)
            T4 = typeof(r)
            new{T1, T2, T3, T4}(CMg0, CFe0, CMn0, Lr, nr, Δr, r, tfinal)
        end
    end
end

@kwdef struct InitialConditions2DMajor{T1, T2, T3, T4, T5} <: InitialConditions
    CMg0::T1
    CFe0::T1
    CMn0::T1
    Lx::T2
    Ly::T2
    nx::T4
    ny::T4
    Δx::T2
    Δy::T2
    x::T5
    y::T5
    grt_position::T3
    grt_boundary::T3
    # grid::NamedTuple{(:x, :y), Tuple{AbstractArray{T2, 2}, AbstractArray{T2, 2}}}
    tfinal::T2
    function InitialConditions2DMajor(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, tfinal::T2, grt_boundary::T3) where {T1 <: AbstractArray{<:Real, 2}, T2 <: Float64, T3 <: Union{AbstractArray{<:Real, 2}, AbstractArray{<:Bool, 2}}}
        if Lx <= 0 || Ly <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) ≠ size(CFe0, 1) || size(CMg0, 1) ≠ size(CMn0, 1) || size(CMg0, 2) ≠ size(CFe0, 2) || size(CMg0, 2) ≠ size(CMn0, 2)
            error("Initial conditions should have the same size.")
        # check that grt_boundary is of the same size as CMg0, CFe0 and CMn0
        elseif size(CMg0, 1) ≠ size(grt_boundary, 1) || size(CMg0, 2) ≠ size(grt_boundary, 2)
            error("grt_boundary should have the same size as CMg0, CFe0 and CMn0.")
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
            grt_position = similar(grt_boundary)
            # create a binary matrix with 1 where grt is present. 0 otherwise. This depends on if CMg0, CFe0 and CMn0 are equal to 0 or not
            grt_position .= (CMg0 .≠ 0) .& (CFe0 .≠ 0) .& (CMn0 .≠ 0)

            T4 = typeof(nx)
            T5 = typeof(x)

            new{T1, T2, T3, T4, T5}(CMg0, CFe0, CMn0, Lx, Ly, nx, ny, Δx, Δy, x, y, grt_position, grt_boundary, tfinal)
        end
    end
end


@kwdef struct InitialConditions3DMajor{T1, T2, T3, T4, T5} <: InitialConditions
    CMg0::T1
    CFe0::T1
    CMn0::T1
    nx::T4
    ny::T4
    nz::T4
    Lx::T2
    Ly::T2
    Lz::T2
    Δx::T2
    Δy::T2
    Δz::T2
    x::T5
    y::T5
    z::T5
    grt_position::T3
    grt_boundary::T3
    tfinal::T2
    function InitialConditions3DMajor(CMg0::T1, CFe0::T1, CMn0::T1, Lx::T2, Ly::T2, Lz::T2, tfinal::T2, grt_boundary::T3) where {T1 <: AbstractArray{<:Real, 3}, T2 <: Float64, T3 <: Union{AbstractArray{<:Real, 3}, AbstractArray{<:Bool, 3}}}
        if Lx <= 0 || Ly <= 0 || Lz <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        elseif size(CMg0, 1) ≠ size(CFe0, 1) || size(CMg0, 1) ≠ size(CMn0, 1) || size(CMg0, 2) ≠ size(CFe0, 2) || size(CMg0, 2) ≠ size(CMn0, 2) || size(CMg0, 3) ≠ size(CFe0, 3) || size(CMg0, 3) ≠ size(CMn0, 3)
            error("Initial conditions should have the same size.")
        elseif size(CMg0, 1) ≠ size(grt_boundary, 1) || size(CMg0, 2) ≠ size(grt_boundary, 2) || size(CMg0, 3) ≠ size(grt_boundary, 3)
            error("grt_boundary should have the same size as CMg0, CFe0 and CMn0.")
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

            # define when grt is present
            grt_position = similar(grt_boundary)
            # create a binary matrix with 1 where grt is present. 0 otherwise. This depends on if CMg0, CFe0 and CMn0 are equal to 0 or not
            grt_position .= (CMg0 .≠ 0) .& (CFe0 .≠ 0) .& (CMn0 .≠ 0)

            T4 = typeof(nx)
            T5 = typeof(x)

            new{T1, T2, T3, T4, T5}(CMg0, CFe0, CMn0, nx, ny, nz, Lx, Ly, Lz, Δx, Δy, Δz, x, y, z, grt_position, grt_boundary, tfinal)
        end
    end
end

"""
    IC1DMajor(CMg0::Array{<:Real}, CFe0::Array{<:Real}, CMn0::Array{<:Real}, Lx::Unitful.Length, x::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(CMg0, 1), stop=Lx), tfinal::Unitful.Time)

Return a structure containing the initial conditions for a 1D diffusion problem. CMg0, CFe0 and CMn0 need to be in molar fraction. Convert `Lx`, `x`, and `tfinal` to µm, µm and Myr, respectively.
"""
function IC1DMajor(;
                   CMg0::AbstractArray{<:Real},
                   CFe0::AbstractArray{<:Real},
                   CMn0::AbstractArray{<:Real},
                   Lx::Unitful.Length,
                   x::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(CMg0, 1), stop=Lx),
                   tfinal::Unitful.Time
                   )

    InitialConditions1DMajor(CMg0, CFe0, CMn0, convert(Float64, ustrip(u"µm", Lx)), ustrip.(u"µm", x), convert(Float64, ustrip(u"Myr",tfinal)))
end


"""
    ICSphMajor(CMg0::AbstractArray{<:Real, 1}, CFe0::AbstractArray{<:Real, 1}, CMn0::AbstractArray{<:Real, 1}, Lr::Unitful.Length, tfinal::Unitful.Time)

Return a structure containing the initial conditions for a spherical diffusion problem. CMg0, CFe0 and CMn0 need to be in molar fraction. Convert `Lr` and `tfinal` to µm and Myr respectively.
"""
function ICSphMajor(;
                    CMg0::AbstractArray{<:Real, 1},
                    CFe0::AbstractArray{<:Real, 1},
                    CMn0::AbstractArray{<:Real, 1},
                    Lr::Unitful.Length,
                    r::AbstractArray{<:Unitful.Length}=range(0u"µm", length=size(CMg0, 1), stop=Lr),
                    tfinal::Unitful.Time
                    )

    InitialConditionsSpherical(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lr)), ustrip.(u"µm", r), convert(Float64,ustrip(u"Myr",tfinal)))
end

"""
    IC2DMajor(;CMg0::AbstractArray{<:Real, 2}, CFe0::AbstractArray{<:Real, 2}, CMn0::AbstractArray{<:Real, 2}, Lx::Unitful.Length, Ly::Unitful.Length, tfinal::Unitful.Time, grt_boundary::AbstractArray{<:Real, 2}=zeros(eltype(CMg0), size(CMg0)...))

Return a structure containing the initial conditions for a 2D diffusion problem. CMg0, CFe0 and CMn0 need to be in molar fraction. Convert `Lx`, `Ly` and `tfinal` to µm, µm and Myr respectively.
grt_boundary is a matrix of the same size as CMg0, CFe0 and CMn0. It is a binary matrix with 1 where the contour of the garnet is present and 0 otherwise. It is used to apply the Dirichlet boundary condition in the model. If not provided, it is equal to zero everywhere.
"""
function IC2DMajor(;
                   CMg0::AbstractArray{<:Real, 2},
                   CFe0::AbstractArray{<:Real, 2},
                   CMn0::AbstractArray{<:Real, 2},
                   Lx::Unitful.Length,
                   Ly::Unitful.Length,
                   tfinal::Unitful.Time,
                   grt_boundary::Union{AbstractArray{<:Real, 2}, AbstractArray{<:Bool, 2}}=zeros(Bool, size(CMg0)...)
                   )

    InitialConditions2DMajor(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"µm", Ly)), convert(Float64,ustrip(u"Myr", tfinal)), grt_boundary)
end

"""
    IC3DMajor(;CMg0::AbstractArray{<:Real, 3}, CFe0::AbstractArray{<:Real, 3}, CMn0::AbstractArray{<:Real, 3}, Lx::Unitful.Length, Ly::Unitful.Length, Lz::Unitful.Length, tfinal::Unitful.Time, grt_boundary::Union{AbstractArray{<:Real, 3}, AbstractArray{<:Bool, 3}}=zeros(Bool, size(CMg0)...))

Return a structure containing the initial conditions for a 3D diffusion problem. CMg0, CFe0 and CMn0 need to be in molar fraction. Convert `Lx`, `Ly`, `Lz` and `tfinal` to µm, µm, µm and Myr respectively.
grt_boundary is a matrix of the same size as CMg0, CFe0 and CMn0. It is a binary matrix with 1 where the contour of the garnet is present and 0 otherwise. It is used to apply the Dirichlet boundary condition in the model. If not provided, it is equal to zero everywhere.
"""
function IC3DMajor(;
                   CMg0::AbstractArray{<:Real, 3},
                   CFe0::AbstractArray{<:Real, 3},
                   CMn0::AbstractArray{<:Real, 3},
                   Lx::Unitful.Length,
                   Ly::Unitful.Length,
                   Lz::Unitful.Length,
                   tfinal::Unitful.Time,
                   grt_boundary::Union{AbstractArray{<:Real, 3}, AbstractArray{<:Bool, 3}}=zeros(Bool, size(CMg0)...)
                   )

    InitialConditions3DMajor(CMg0, CFe0, CMn0, convert(Float64,ustrip(u"µm", Lx)), convert(Float64,ustrip(u"µm", Ly)), convert(Float64,ustrip(u"µm", Lz)), convert(Float64,ustrip(u"Myr", tfinal)), grt_boundary)
end


@kwdef struct Domain1DMajor{T1, T2, T3, T4, T5, T_tuplediffdata} <: Domain
    IC::T4
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    diffcoef::Int
    D0_data::T_tuplediffdata  # tuple of diffusion coefficients data
    D0::Matrix{T2}  # diffusion coefficients
    D::NamedTuple{(:DMgMg, :DMgFe, :DMgMn, :DFeMg, :DFeFe, :DFeMn, :DMnMg, :DMnFe, :DMnMn),
                  NTuple{9, Vector{T2}}}  # matrix of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T5
    u0::Matrix{T2}
    tfinal_ad::T2
    time_update_ad::T1
    bc_neumann::T3
    function Domain1DMajor(IC::InitialConditions1DMajor, T::T1, P::T1, time_update::T1, fugacity_O2::T1, diffcoef::Int, bc_neumann::T3) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}, T3 <: Tuple}

        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, Δx, tfinal, Lx, CMg0, CFe0, CMn0 = IC

        # define the trace diffusion coefficients based on the chosen dataset
        if diffcoef == 1
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chakraborty1992)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chakraborty1992)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chakraborty1992)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn)
        elseif diffcoef == 2
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Carlson2006)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Carlson2006)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Carlson2006)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Carlson2006)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        elseif diffcoef == 3
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chu2015)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chu2015)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chu2015)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Chu2015)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        end

        T_tuplediffdata = typeof(D0_data)

        T2 = eltype(CMg0)
        D0 = similar(CMg0, (4, nx))

        # iterate through D0 and compute the initial diffusion coefficients
        for j in axes(D0, 2)
            D0_view = @view D0[:, j]

            T_K = (T[1] + 273.15) * u"K"
            P_kbar = P[1] * u"kbar"
            fO2 = (fugacity_O2[1])NoUnits

            # compute the diffusion coefficients for each point
            D_update!(D0_view, T_K, P_kbar, diffcoef, CMg0[j], CFe0[j], CMn0[j], D0_data, fO2)  # compute initial diffusion coefficients
        end

        D = (DMgMg = similar(CMg0, nx),
             DMgFe = similar(CMg0, nx),
             DMgMn = similar(CMg0, nx),
             DFeMg = similar(CMg0, nx),
             DFeFe = similar(CMg0, nx),
             DFeMn = similar(CMg0, nx),
             DMnMg = similar(CMg0, nx),
             DMnFe = similar(CMg0, nx),
             DMnMn = similar(CMg0, nx))  # matrix of interdiffusion coefficients

        for i in eachindex(D)
            D[i] .= 0  # copy the initial diffusion coefficients to the interdiffusion coefficients
        end

        u0 = similar(CMg0, (nx, 3))
        u0[:,1] .= CMg0
        u0[:,2] .= CFe0
        u0[:,3] .= CMn0

        # # find index of the first column of non-zero values in D0
        # first_nonzero_col = findfirst(!iszero, D0[1, :])
        # D0_mean = mean(D0[:, first_nonzero_col])  # mean diffusion coefficient

        # L_charact = copy(Lx)  # characteristic length
        # D_charact = D0_mean  # characteristic
        # t_charact = L_charact^2 / D0_mean  # characteristic time

        # Now fix t to be 1 and make D_charact = L_charact^2 / t_charact
        L_charact = copy(Lx)  # characteristic length
        t_charact = 1.0
        D_charact = L_charact^2 / t_charact

        Δxad_ = 1 ./ (Δx ./ L_charact)  # inverse of nondimensionalised Δx
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T4 = typeof(IC)
        T5 = typeof(Δxad_)

        new{T1, T2, T3, T4, T5, T_tuplediffdata}(IC, T, P, fugacity_O2, time_update, diffcoef, D0_data, D0, D, L_charact, D_charact, t_charact, Δxad_, u0, tfinal_ad, time_update_ad, bc_neumann)
    end
end

@kwdef struct DomainSphericalMajor{T1, T2, T3, T_tuplediffdata} <: Domain
    IC::T3
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    diffcoef::Int
    D0_data::T_tuplediffdata  # tuple of diffusion coefficients data
    D0::Matrix{T2}
    D::NamedTuple{(:DMgMg, :DMgFe, :DMgMn, :DFeMg, :DFeFe, :DFeMn, :DMnMg, :DMnFe, :DMnMn),
                  NTuple{9, Vector{T2}}}  # matrix of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δr_ad::Vector{T2}
    Δr_ad_::Vector{T2}
    r_ad::Vector{T2}
    u0::Matrix{T2}
    tfinal_ad::T2
    time_update_ad::T1
    function DomainSphericalMajor(IC::InitialConditionsSpherical, T::T1, P::T1, time_update::T1, fugacity_O2::T1, diffcoef::Int) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}

        #check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nr, Δr, r, tfinal, Lr, CMg0, CFe0, CMn0 = IC

        # define the trace diffusion coefficients based on the chosen dataset
        if diffcoef == 1
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chakraborty1992)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chakraborty1992)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chakraborty1992)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn)
        elseif diffcoef == 2
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Carlson2006)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Carlson2006)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Carlson2006)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Carlson2006)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        elseif diffcoef == 3
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chu2015)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chu2015)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chu2015)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Chu2015)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        end

        T_tuplediffdata = typeof(D0_data)

        T2 = eltype(CMg0)

        D0 = similar(CMg0, (4, nr))

        # iterate through D0 and compute the initial diffusion coefficients
        for j in axes(D0, 2)
            D0_view = @view D0[:, j]

            T_K = (T[1] + 273.15) * u"K"
            P_kbar = P[1] * u"kbar"
            fO2 = (fugacity_O2[1])NoUnits

            # compute the diffusion coefficients for each point
            D_update!(D0_view, T_K, P_kbar, diffcoef, CMg0[j], CFe0[j], CMn0[j], D0_data, fO2)  # compute initial diffusion coefficients
        end

        D = (DMgMg = similar(CMg0, nr),
             DMgFe = similar(CMg0, nr),
             DMgMn = similar(CMg0, nr),
             DFeMg = similar(CMg0, nr),
             DFeFe = similar(CMg0, nr),
             DFeMn = similar(CMg0, nr),
             DMnMg = similar(CMg0, nr),
             DMnFe = similar(CMg0, nr),
             DMnMn = similar(CMg0, nr))  # matrix of interdiffusion coefficients

        for i in eachindex(D)
            D[i] .= 0
        end

        u0::Matrix{T2} = similar(CMg0, (nr, 3))
        u0[:,1] .= CMg0
        u0[:,2] .= CFe0
        u0[:,3] .= CMn0

        # find index of the first column of non-zero values in D0
        # first_nonzero_col = findfirst(!iszero, D0[1, :])
        # D0_mean = mean(D0[:, first_nonzero_col])  # mean diffusion coefficient

        # L_charact = copy(Lr)  # characteristic length
        # D_charact = D0_mean  # characteristic
        # t_charact = L_charact^2 / D_charact  # characteristic time

        # Now fix t to be 1 and make D_charact = L_charact^2 / t_charact
        L_charact = Lr  # characteristic length
        t_charact = 1.0
        D_charact = L_charact^2 / t_charact

        Δr_ad = Δr ./ L_charact  # nondimensionalised Δr
        Δr_ad_ = 1 ./ Δr_ad  # inverse of nondimensionalised Δr
        r_ad = r ./ L_charact  # nondimensionalised radius
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T3 = typeof(IC)

        new{T1, T2, T3, T_tuplediffdata}(IC, T, P, fugacity_O2, time_update, diffcoef, D0_data, D0, D, L_charact, D_charact, t_charact, Δr_ad, Δr_ad_, r_ad, u0, tfinal_ad, time_update_ad)
    end
end

@kwdef struct Domain2DMajor{T1, T2, T3, T4, T5, T6, T_tuplediffdata} <: Domain
    IC::T6
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    diffcoef::Int
    D0_data::T_tuplediffdata  # tuple of diffusion coefficients data
    D0::T5
    D::T4  # matrix of interdiffusion coefficients
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    Δyad_::T2
    u0::T5
    tfinal_ad::T2
    time_update_ad::T1
    function Domain2DMajor(IC::InitialConditions2DMajor, T::T1, P::T1, time_update::T1, fugacity_O2::T1, diffcoef::Int) where {T1 <: Union{Float64, AbstractArray{Float64, 1}}}

        # check that T, P and time_update have the same size
        if size(T, 1) ≠ size(P, 1) || size(T, 1) ≠ size(time_update, 1)
            error("T, P and time_update should have the same size.")
        end

        @unpack nx, ny, Δx, Δy, tfinal, Lx, CMg0, CFe0, CMn0, grt_position, grt_boundary = IC

        # define the trace diffusion coefficients based on the chosen dataset
        if diffcoef == 1
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chakraborty1992)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chakraborty1992)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chakraborty1992)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Fe) # fake Grt_Ca to make parallelStencil happy!!
        elseif diffcoef == 2
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Carlson2006)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Carlson2006)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Carlson2006)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Carlson2006)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        elseif diffcoef == 3
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chu2015)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chu2015)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chu2015)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Chu2015)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        end

        T_tuplediffdata = typeof(D0_data)

        D0 = similar(CMg0, (4, nx, ny))

        T_K = (T[1] + 273.15) * u"K"
        P_kbar = P[1] * u"kbar"
        fO2 = (fugacity_O2[1])NoUnits

        @parallel D_update_2D!(D0, T_K, P_kbar, diffcoef, CMg0, CFe0, CMn0, D0_data, fO2, grt_position, grt_boundary)

        D = (DMgMg = similar(CMg0, (nx, ny)),
             DMgFe = similar(CMg0, (nx, ny)),
             DMgMn = similar(CMg0, (nx, ny)),
             DFeMg = similar(CMg0, (nx, ny)),
             DFeFe = similar(CMg0, (nx, ny)),
             DFeMn = similar(CMg0, (nx, ny)),
             DMnMg = similar(CMg0, (nx, ny)),
             DMnFe = similar(CMg0, (nx, ny)),
             DMnMn = similar(CMg0, (nx, ny)))  # matrix of interdiffusion coefficients

        for i in eachindex(D)
            D[i] .= 0
        end

        u0 = similar(CMg0, (nx, ny, 3))
        u0[:, :, 1] .= CMg0
        u0[:, :, 2] .= CFe0
        u0[:, :, 3] .= CMn0

        # Now fix t to be 1 and make D_charact = L_charact^2 / t_charact
        L_charact = Lx  # characteristic length
        t_charact = 1.0
        D_charact = L_charact^2 / t_charact

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        Δyad_ = 1 / (Δy / L_charact)  # inverse of nondimensionalised Δy
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T2 = typeof(t_charact)
        T3 = typeof(D0)
        T4 = typeof(D)
        T5 = typeof(u0)
        T6 = typeof(IC)

        new{T1, T2, T3, T4, T5, T6, T_tuplediffdata}(IC, T, P, fugacity_O2, time_update, diffcoef, D0_data, D0, D, L_charact, D_charact, t_charact, Δxad_, Δyad_, u0, tfinal_ad, time_update_ad)
    end
end

@kwdef struct Domain3DMajor{T1, T2, T3, T4, T5, T6, T_tuplediffdata} <: Domain
    IC::T6
    T::T1
    P::T1
    fugacity_O2::T1
    time_update::T1
    diffcoef::Int
    D0_data::T_tuplediffdata  # tuple of diffusion coefficients data
    D0::T3
    D::T4
    L_charact::T2
    D_charact::T2
    t_charact::T2
    Δxad_::T2
    Δyad_::T2
    Δzad_::T2
    u0::T5
    tfinal_ad::T2
    time_update_ad::T1
    function Domain3DMajor(IC::InitialConditions3DMajor, T::T1, P::T1, time_update::T1, fugacity_O2::T1, diffcoef::Int) where {T1 <: Union{Float64, Array{Float64, 1}}}

        @unpack nx, ny, nz, Δx, Δy, Δz, tfinal, Lx, CMg0, CFe0, CMn0, grt_position, grt_boundary = IC

        # define the trace diffusion coefficients based on the chosen dataset
        if diffcoef == 1
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chakraborty1992)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chakraborty1992)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chakraborty1992)
            Grt_Ca = SetChemicalDiffusion(Grt_Fe_Chakraborty1992) # fake one to make GPU happy

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        elseif diffcoef == 2
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Carlson2006)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Carlson2006)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Carlson2006)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Carlson2006)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        elseif diffcoef == 3
            Grt_Mg = SetChemicalDiffusion(Grt_Mg_Chu2015)
            Grt_Fe = SetChemicalDiffusion(Grt_Fe_Chu2015)
            Grt_Mn = SetChemicalDiffusion(Grt_Mn_Chu2015)
            Grt_Ca = SetChemicalDiffusion(Grt_Ca_Chu2015)

            D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn, Grt_Ca=Grt_Ca)
        end

        T_tuplediffdata = typeof(D0_data)

        D0 = similar(CMg0, (4, nx, ny, nz))

        T_K = (T[1] + 273.15) * u"K"
        P_kbar = P[1] * u"kbar"
        fO2 = (fugacity_O2[1])NoUnits

        @parallel D_update_3D!(D0, T_K, P_kbar, diffcoef, CMg0, CFe0, CMn0, D0_data, fO2, grt_position, grt_boundary)

        D = (DMgMg = similar(CMg0, (nx, ny, nz)),
             DMgFe = similar(CMg0, (nx, ny, nz)),
             DMgMn = similar(CMg0, (nx, ny, nz)),
             DFeMg = similar(CMg0, (nx, ny, nz)),
             DFeFe = similar(CMg0, (nx, ny, nz)),
             DFeMn = similar(CMg0, (nx, ny, nz)),
             DMnMg = similar(CMg0, (nx, ny, nz)),
             DMnFe = similar(CMg0, (nx, ny, nz)),
             DMnMn = similar(CMg0, (nx, ny, nz)))  # matrix of interdiffusion coefficients

        for i in eachindex(D)
            D[i] .= 0
        end

        u0 = similar(CMg0, (nx, ny, nz, 3))
        u0[:, :, :, 1] .= CMg0
        u0[:, :, :, 2] .= CFe0
        u0[:, :, :, 3] .= CMn0

        # Now fix t to be 1 and make D_charact = L_charact^2 / t_charact
        L_charact = Lx  # characteristic length
        t_charact = 1.0
        D_charact = L_charact^2 / t_charact

        Δxad_ = 1 / (Δx / L_charact)  # inverse of nondimensionalised Δx
        Δyad_ = 1 / (Δy / L_charact)  # inverse of nondimensionalised Δy
        Δzad_ = 1 / (Δz / L_charact)  # inverse of nondimensionalised Δz
        tfinal_ad = tfinal / t_charact  # nondimensionalised total time
        time_update_ad = time_update ./ t_charact  # nondimensionalised time update

        T2 = typeof(t_charact)
        T3 = typeof(D0)
        T4 = typeof(D)
        T5 = typeof(u0)
        T6 = typeof(IC)

        new{T1, T2, T3, T4, T5, T6, T_tuplediffdata}(IC, T, P, fugacity_O2, time_update, diffcoef, D0_data, D0, D, L_charact, D_charact, t_charact, Δxad_, Δyad_, Δzad_, u0, tfinal_ad, time_update_ad)
    end
end


"""
    Domain(IC::InitialConditions1D, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr"; bc_neumann::Tuple=(false, false))

When applied to 1D initial conditions, define corresponding 1D domain. `bc_neumann` can be used to define Neumann boundary conditions on the left or right side of the domain if set to true.
"""
function Domain(IC::InitialConditions1DMajor, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr", fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa";diffcoef::Symbol=:CG92, bc_neumann::Tuple=(false, false))

    if diffcoef == :CG92
        diffcoef = 1
    elseif diffcoef == :C06
        diffcoef = 2
    elseif diffcoef == :CA15
        diffcoef = 3
    else
        error("Unknown diffusion coefficient. Use :CG92 for Chakraborty and Ganguly (1992), :C06 for Carlson 2006 and :CA15 for Chu and Ague 2015.")
    end

    Domain1DMajor(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), Int(diffcoef), bc_neumann)
end

"""
    Domain(IC::InitialConditionsSpherical, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to spherical initial conditions, define corresponding spherical domain. Assume that the center of the grain is on the left side.
"""
function Domain(IC::InitialConditionsSpherical, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr", fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa"; diffcoef::Symbol=:CG92)

    if diffcoef == :CG92
        diffcoef = 1
    elseif diffcoef == :C06
        diffcoef = 2
    elseif diffcoef == :CA15
        diffcoef = 3
    else
        error("Unknown diffusion coefficient. Use :CG92 for Chakraborty and Ganguly (1992), :C06 for Carlson 2006 and :CA15 for Chu and Ague 2015.")
    end


    DomainSphericalMajor(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), Int(diffcoef))
end

"""
    Domain(IC::IC2DMajor, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 2D initial conditions, define corresponding 2D domain.
"""
function Domain(IC::InitialConditions2DMajor, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr", fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa"; diffcoef::Symbol=:CG92)

    if diffcoef == :CG92
        diffcoef = 1
    elseif diffcoef == :C06
        diffcoef = 2
    elseif diffcoef == :CA15
        diffcoef = 3
    else
        error("Unknown diffusion coefficient. Use :CG92 for Chakraborty and Ganguly (1992), :C06 for Carlson 2006 and :CA15 for Chu and Ague 2015.")
    end

    Domain2DMajor(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), Int(diffcoef))
end

"""
    Domain(IC::IC3DMajor, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr")

When applied to 3D initial conditions, define corresponding 3D domain.
"""
function Domain(IC::InitialConditions3DMajor, T::Union{Unitful.Temperature,Array{<:Unitful.Temperature{<:Real}, 1}}, P::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}, time_update::Union{Unitful.Time,Array{<:Unitful.Time{<:Real}, 1}}=0u"Myr", fugacity_O2::Union{Unitful.Pressure,Array{<:Unitful.Pressure{<:Real}, 1}}=ones(size(P)) .* 1e-25u"Pa"; diffcoef::Symbol=:CG92)

    if diffcoef == :CG92
        diffcoef = 1
    elseif diffcoef == :C06
        diffcoef = 2
    elseif diffcoef == :CA15
        diffcoef = 3
    else
        error("Unknown diffusion coefficient. Use :CG92 for Chakraborty and Ganguly (1992), :C06 for Carlson 2006 and :CA15 for Chu and Ague 2015.")
    end

    Domain3DMajor(IC, convert.(Float64,ustrip.(u"°C", T)), convert.(Float64,ustrip.(u"kbar", P)), convert.(Float64,ustrip.(u"Myr", time_update)), convert.(Float64,ustrip.(u"Pa", fugacity_O2)), Int(diffcoef))
end