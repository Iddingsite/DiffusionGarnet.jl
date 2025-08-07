import Base.@propagate_inbounds

function Diffusion_coef_1D_major!(D, CMg, CFe, CMn, D0, D_charact, domain, t)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @unpack diffcoef, D0_data, T, P, fugacity_O2, time_update_ad = domain

    D_charact_ = 1 / D_charact

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, I) = D0[1, I] * CMg[I] + D0[2, I] * CFe[I] + D0[3, I] * CMn[I] +
    D0[4, I] * (1 - CMg[I] - CFe[I] - CMn[I])


    # end-member unit-cell dimensions for C12 and CA15
    a0_Fe = 1.1525
    a0_Mg = 1.1456
    a0_Mn = 1.1614
    a0_Ca = 1.1852

    # find current T, P, and fugacity_O2
    # to do this, we need to find the lower bound of t in time_update_ad
    index = findfirst(x -> x >= t, time_update_ad)
    if index === nothing
        index = length(time_update_ad)
    end

    P_kbar = P[index] * 1u"kbar"
    T_K = (T[index] + 273.15) * 1u"K"
    fO2 = (fugacity_O2[index])NoUnits

    # @inbounds for I in eachindex(DMgMg)
    for I in eachindex(DMgMg)

        # there is a composition dependence in the self-diffusion coefficients for C12 and CA15
        if diffcoef == 2 || diffcoef == 3

            X = (CFe[I] * a0_Fe + CMg[I] * a0_Mg + CMn[I] * a0_Mn + (1 - (CMg[I] + CFe[I] + CMn[I])) * a0_Ca)NoUnits

            D0[1, I] = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
            D0[2, I] = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
            D0[3, I] = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
            D0[4, I] = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
        end

        sum_D_ = 1 / (sum_D(CMg, CFe, CMn, D0, I))

        DMgMg[I] = (D0[1, I] - D0[1, I] * CMg[I] * sum_D_ * (D0[1, I] - D0[end, I])) * D_charact_
        DMgFe[I] = (      - D0[1, I] * CMg[I] * sum_D_ * (D0[2, I] - D0[end, I])) * D_charact_
        DMgMn[I] = (      - D0[1, I] * CMg[I] * sum_D_ * (D0[3, I] - D0[end, I])) * D_charact_

        DFeMg[I] = (      - D0[2, I] * CFe[I] * sum_D_ * (D0[1, I] - D0[end, I])) * D_charact_
        DFeFe[I] = (D0[2, I] - D0[2, I] * CFe[I] * sum_D_ * (D0[2, I] - D0[end, I])) * D_charact_
        DFeMn[I] = (      - D0[2, I] * CFe[I] * sum_D_ * (D0[3, I] - D0[end, I])) * D_charact_

        DMnMg[I] = (      - D0[3, I] * CMn[I] * sum_D_ * (D0[1, I] - D0[end, I])) * D_charact_
        DMnFe[I] = (      - D0[3, I] * CMn[I] * sum_D_ * (D0[2, I] - D0[end, I])) * D_charact_
        DMnMn[I] = (D0[3, I] - D0[3, I] * CMn[I] * sum_D_ * (D0[3, I] - D0[end, I])) * D_charact_
    end
end

function stencil_diffusion_1D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_, bc_neumann)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δxad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * (Δxad_[ix])
    @propagate_inbounds @inline neumann_left(D, C, ix, Δxad_) = - D[ix] * (C[ix]-C[ix+1]) * (Δxad_[ix])
    @propagate_inbounds @inline neumann_right(D, C, ix, Δxad_) = D[ix] * (C[ix-1]-C[ix]) * (Δxad_[ix-1])

    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, Δxad_)
        Δxad_centered = (Δxad_[ix] + Δxad_[ix-1]) / 2

        dtC[ix] = (qx(D1,C1,ix,Δxad_) - qx(D1,C1,ix-1,Δxad_)) * Δxad_centered +
                  (qx(D2,C2,ix,Δxad_) - qx(D2,C2,ix-1,Δxad_)) * Δxad_centered +
                  (qx(D3,C3,ix,Δxad_) - qx(D3,C3,ix-1,Δxad_)) * Δxad_centered
    end

    for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, Δxad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, Δxad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, Δxad_)
        end
    end

    # neumann boundary conditions
    if bc_neumann[1] == true
        ix = 1
        Δxad_centered = Δxad_[ix]

        dtCMg[1] = qx(DMgMg,CMg,ix,Δxad_) * Δxad_centered +
                   qx(DMgFe,CFe,ix,Δxad_) * Δxad_centered +
                   qx(DMgMn,CMn,ix,Δxad_) * Δxad_centered
        dtCFe[1] = qx(DFeMg,CMg,ix,Δxad_) * Δxad_centered +
                   qx(DFeFe,CFe,ix,Δxad_) * Δxad_centered +
                   qx(DFeMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMn[1] = qx(DMnMg,CMg,ix,Δxad_) * Δxad_centered +
                   qx(DMnFe,CFe,ix,Δxad_) * Δxad_centered +
                   qx(DMnMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMg[1] += neumann_left(DMgMg,CMg,ix,Δxad_) * Δxad_centered +
                    neumann_left(DMgFe,CFe,ix,Δxad_) * Δxad_centered +
                    neumann_left(DMgMn,CMn,ix,Δxad_) * Δxad_centered
        dtCFe[1] += neumann_left(DFeMg,CMg,ix,Δxad_) * Δxad_centered +
                    neumann_left(DFeFe,CFe,ix,Δxad_) * Δxad_centered +
                    neumann_left(DFeMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMn[1] += neumann_left(DMnMg,CMg,ix,Δxad_) * Δxad_centered +
                    neumann_left(DMnFe,CFe,ix,Δxad_) * Δxad_centered +
                    neumann_left(DMnMn,CMn,ix,Δxad_) * Δxad_centered
    end

    if bc_neumann[2] == true
        ix = last_index(dtCMg)
        Δxad_centered = Δxad_[ix]

        dtCMg[end] = - qx(DMgMg,CMg,ix,Δxad_) * Δxad_centered -
                       qx(DMgFe,CFe,ix,Δxad_) * Δxad_centered -
                       qx(DMgMn,CMn,ix,Δxad_) * Δxad_centered
        dtCFe[end] = - qx(DFeMg,CMg,ix,Δxad_) * Δxad_centered -
                       qx(DFeFe,CFe,ix,Δxad_) * Δxad_centered -
                       qx(DFeMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMn[end] = - qx(DMnMg,CMg,ix,Δxad_) * Δxad_centered -
                       qx(DMnFe,CFe,ix,Δxad_) * Δxad_centered -
                       qx(DMnMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMg[end] += neumann_right(DMgMg,CMg,ix,Δxad_) * Δxad_centered +
                      neumann_right(DMgFe,CFe,ix,Δxad_) * Δxad_centered +
                      neumann_right(DMgMn,CMn,ix,Δxad_) * Δxad_centered
        dtCFe[end] += neumann_right(DFeMg,CMg,ix,Δxad_) * Δxad_centered +
                      neumann_right(DFeFe,CFe,ix,Δxad_) * Δxad_centered +
                      neumann_right(DFeMn,CMn,ix,Δxad_) * Δxad_centered
        dtCMn[end] += neumann_right(DMnMg,CMg,ix,Δxad_) * Δxad_centered +
                      neumann_right(DMnFe,CFe,ix,Δxad_) * Δxad_centered +
                      neumann_right(DMnMn,CMn,ix,Δxad_) * Δxad_centered
    end
end


function semi_discretisation_diffusion_cartesian(du::T,u::T,p,t) where T <: AbstractArray{<:Real, 2}

    @unpack D, D0, D_charact, Δxad_, bc_neumann = p.domain

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    Diffusion_coef_1D_major!(D, CMg, CFe ,CMn, D0, D_charact, p.domain, t)

    # semi-discretization
    stencil_diffusion_1D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_, bc_neumann)
end

