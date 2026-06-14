import Base.@propagate_inbounds

@parallel_indices (ix, iy, iz) function Diffusion_coef_3D_major!(D, CMg, CFe ,CMn, D0, D_charact_, grt_position, diffcoef, D0_data, T_K, P_kbar, fO2)

    @muladd @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0Mg, D0Fe, D0Mn, D0Ca, ix, iy, iz) = D0Mg * CMg[ix, iy, iz] + D0Fe * CFe[ix, iy, iz] + D0Mn * CMn[ix, iy, iz] +
        D0Ca * (1 - CMg[ix, iy, iz] - CFe[ix, iy, iz] - CMn[ix, iy, iz])

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @inbounds begin
        if ix>1 && ix<size(DMgMg,1) && iy>1 && iy<size(DMgMg,2) && iz>1 && iz<size(DMgMg,3)
            if isone(grt_position[ix,iy,iz])

                D0Mg = D0[1]
                D0Fe = D0[2]
                D0Mn = D0[3]
                D0Ca = D0[4]

                if diffcoef == 2 || diffcoef == 3

                    # end-member unit-cell dimensions for C12 and CA15
                    a0_Fe = 1.1525
                    a0_Mg = 1.1456
                    a0_Mn = 1.1614
                    a0_Ca = 1.1852

                    X = CFe[ix,iy,iz] * a0_Fe + CMg[ix,iy,iz] * a0_Mg + CMn[ix,iy,iz] * a0_Mn + (1 - (CMg[ix,iy,iz] + CFe[ix,iy,iz] + CMn[ix,iy,iz])) * a0_Ca
                    X = convert(eltype(CMg), X)NoUnits

                    D0Mg = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fO2, X = X))))
                    D0Fe = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fO2, X = X))))
                    D0Mn = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fO2, X = X))))
                    D0Ca = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = fO2, X = X))))
                end

                # prevent reading 3 times the same value on GPU
                cMg = CMg[ix,iy,iz]
                cFe = CFe[ix,iy,iz]
                cMn = CMn[ix,iy,iz]

                # factor that repeats 3 times for each diffusion coefficient
                sum_D_ = inv(sum_D(CMg,CFe,CMn,D0Mg,D0Fe,D0Mn,D0Ca,ix,iy,iz))

                sum_DMg = sum_D_ * (D0Mg - D0Ca)
                sum_DFe = sum_D_ * (D0Fe - D0Ca)
                sum_DMn = sum_D_ * (D0Mn - D0Ca)

                # factor that repeats 3 times for each diffusion coefficient
                D0MgCMg = D0Mg * cMg
                D0FeCFe = D0Fe * cFe
                D0MnCMn = D0Mn * cMn

                DMgMg[ix,iy,iz] = (D0Mg - D0MgCMg * sum_DMg) * D_charact_
                DMgFe[ix,iy,iz] = (     - D0MgCMg * sum_DFe) * D_charact_
                DMgMn[ix,iy,iz] = (     - D0MgCMg * sum_DMn) * D_charact_
                DFeMg[ix,iy,iz] = (     - D0FeCFe * sum_DMg) * D_charact_
                DFeFe[ix,iy,iz] = (D0Fe - D0FeCFe * sum_DFe) * D_charact_
                DFeMn[ix,iy,iz] = (     - D0FeCFe * sum_DMn) * D_charact_
                DMnMg[ix,iy,iz] = (     - D0MnCMn * sum_DMg) * D_charact_
                DMnFe[ix,iy,iz] = (     - D0MnCMn * sum_DFe) * D_charact_
                DMnMn[ix,iy,iz] = (D0Mn - D0MnCMn * sum_DMn) * D_charact_

            end
        end
    end

    return
end


@parallel_indices (ix, iy, iz) function stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, position_Grt, Grt_boundaries, Δxad_², Δyad_², Δzad_²)

    # raw weighted differences — Δ² factor applied once at divergence level
    @propagate_inbounds @inline av_D_x(D, ix, iy, iz) = (D[ix,iy,iz] + D[ix+1,iy,iz]) / 2
    @propagate_inbounds @inline av_D_y(D, ix, iy, iz) = (D[ix,iy,iz] + D[ix,iy+1,iz]) / 2
    @propagate_inbounds @inline av_D_z(D, ix, iy, iz) = (D[ix,iy,iz] + D[ix,iy,iz+1]) / 2
    @propagate_inbounds @inline qx(D, C, ix, iy, iz) = av_D_x(D, ix, iy, iz) * (C[ix+1,iy,iz]-C[ix,iy,iz])
    @propagate_inbounds @inline qy(D, C, ix, iy, iz) = av_D_y(D, ix, iy, iz) * (C[ix,iy+1,iz]-C[ix,iy,iz])
    @propagate_inbounds @inline qz(D, C, ix, iy, iz) = av_D_z(D, ix, iy, iz) * (C[ix,iy,iz+1]-C[ix,iy,iz])


    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, iy, iz, Δxad_², Δyad_², Δzad_²)
        dtC[ix,iy,iz] = (qx(D1,C1,ix,iy,iz) - qx(D1,C1,ix-1,iy,iz)) * Δxad_² +
                        (qx(D2,C2,ix,iy,iz) - qx(D2,C2,ix-1,iy,iz)) * Δxad_² +
                        (qx(D3,C3,ix,iy,iz) - qx(D3,C3,ix-1,iy,iz)) * Δxad_² +
                        (qy(D1,C1,ix,iy,iz) - qy(D1,C1,ix,iy-1,iz)) * Δyad_² +
                        (qy(D2,C2,ix,iy,iz) - qy(D2,C2,ix,iy-1,iz)) * Δyad_² +
                        (qy(D3,C3,ix,iy,iz) - qy(D3,C3,ix,iy-1,iz)) * Δyad_² +
                        (qz(D1,C1,ix,iy,iz) - qz(D1,C1,ix,iy,iz-1)) * Δzad_² +
                        (qz(D2,C2,ix,iy,iz) - qz(D2,C2,ix,iy,iz-1)) * Δzad_² +
                        (qz(D3,C3,ix,iy,iz) - qz(D3,C3,ix,iy,iz-1)) * Δzad_²
    end

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    # iterate inside the arrays
    if ix>1 && ix<size(dtCMg,1) && iy>1 && iy<size(dtCMg,2) && iz>1 && iz<size(dtCMg,3)
        if isone(position_Grt[ix,iy,iz]) && iszero(Grt_boundaries[ix,iy,iz])
            @inbounds update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, iy, iz, Δxad_², Δyad_², Δzad_²)
            @inbounds update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, iy, iz, Δxad_², Δyad_², Δzad_²)
            @inbounds update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, iy, iz, Δxad_², Δyad_², Δzad_²)

            # first order Neumann if inclusions
            # south
            if iszero(position_Grt[ix-1,iy,iz])
                @inbounds dtCMg[ix,iy,iz] -= - qx(DMgMg,CMg,ix-1,iy,iz) * Δxad_² -
                                  qx(DMgFe,CFe,ix-1,iy,iz) * Δxad_² -
                                  qx(DMgMn,CMn,ix-1,iy,iz) * Δxad_²
                @inbounds dtCFe[ix,iy,iz] -= - qx(DFeMg,CMg,ix-1,iy,iz) * Δxad_² -
                                  qx(DFeFe,CFe,ix-1,iy,iz) * Δxad_² -
                                  qx(DFeMn,CMn,ix-1,iy,iz) * Δxad_²
                @inbounds dtCMn[ix,iy,iz] -= - qx(DMnMg,CMg,ix-1,iy,iz) * Δxad_² -
                                  qx(DMnFe,CFe,ix-1,iy,iz) * Δxad_² -
                                  qx(DMnMn,CMn,ix-1,iy,iz) * Δxad_²
            end
            # north
            if iszero(position_Grt[ix+1,iy,iz])
                @inbounds dtCMg[ix,iy,iz] -= qx(DMgMg,CMg,ix,iy,iz) * Δxad_² +
                                qx(DMgFe,CFe,ix,iy,iz) * Δxad_² +
                                qx(DMgMn,CMn,ix,iy,iz) * Δxad_²
                @inbounds dtCFe[ix,iy,iz] -= qx(DFeMg,CMg,ix,iy,iz) * Δxad_² +
                                qx(DFeFe,CFe,ix,iy,iz) * Δxad_² +
                                qx(DFeMn,CMn,ix,iy,iz) * Δxad_²
                @inbounds dtCMn[ix,iy,iz] -= qx(DMnMg,CMg,ix,iy,iz) * Δxad_² +
                                qx(DMnFe,CFe,ix,iy,iz) * Δxad_² +
                                qx(DMnMn,CMn,ix,iy,iz) * Δxad_²
            end
            # west
            if iszero(position_Grt[ix,iy-1,iz])
                @inbounds dtCMg[ix,iy,iz] -= - qy(DMgMg,CMg,ix,iy-1,iz) * Δyad_² -
                                  qy(DMgFe,CFe,ix,iy-1,iz) * Δyad_² -
                                  qy(DMgMn,CMn,ix,iy-1,iz) * Δyad_²
                @inbounds dtCFe[ix,iy,iz] -= - qy(DFeMg,CMg,ix,iy-1,iz) * Δyad_² -
                                  qy(DFeFe,CFe,ix,iy-1,iz) * Δyad_² -
                                  qy(DFeMn,CMn,ix,iy-1,iz) * Δyad_²
                @inbounds dtCMn[ix,iy,iz] -= - qy(DMnMg,CMg,ix,iy-1,iz) * Δyad_² -
                                  qy(DMnFe,CFe,ix,iy-1,iz) * Δyad_² -
                                  qy(DMnMn,CMn,ix,iy-1,iz) * Δyad_²
            end
            # east
            if iszero(position_Grt[ix,iy+1,iz])
                @inbounds dtCMg[ix,iy,iz] -= qy(DMgMg,CMg,ix,iy,iz) * Δyad_² +
                                qy(DMgFe,CFe,ix,iy,iz) * Δyad_² +
                                qy(DMgMn,CMn,ix,iy,iz) * Δyad_²
                @inbounds dtCFe[ix,iy,iz] -= qy(DFeMg,CMg,ix,iy,iz) * Δyad_² +
                                qy(DFeFe,CFe,ix,iy,iz) * Δyad_² +
                                qy(DFeMn,CMn,ix,iy,iz) * Δyad_²
                @inbounds dtCMn[ix,iy,iz] -= qy(DMnMg,CMg,ix,iy,iz) * Δyad_² +
                                qy(DMnFe,CFe,ix,iy,iz) * Δyad_² +
                                qy(DMnMn,CMn,ix,iy,iz) * Δyad_²
            end
            # bottom
            if iszero(position_Grt[ix,iy,iz-1])
                @inbounds dtCMg[ix,iy,iz] -= - qz(DMgMg,CMg,ix,iy,iz-1) * Δzad_² -
                                  qz(DMgFe,CFe,ix,iy,iz-1) * Δzad_² -
                                  qz(DMgMn,CMn,ix,iy,iz-1) * Δzad_²
                @inbounds dtCFe[ix,iy,iz] -= - qz(DFeMg,CMg,ix,iy,iz-1) * Δzad_² -
                                  qz(DFeFe,CFe,ix,iy,iz-1) * Δzad_² -
                                  qz(DFeMn,CMn,ix,iy,iz-1) * Δzad_²
                @inbounds dtCMn[ix,iy,iz] -= - qz(DMnMg,CMg,ix,iy,iz-1) * Δzad_² -
                                  qz(DMnFe,CFe,ix,iy,iz-1) * Δzad_² -
                                  qz(DMnMn,CMn,ix,iy,iz-1) * Δzad_²
            end
            # top
            if iszero(position_Grt[ix,iy,iz+1])
                @inbounds dtCMg[ix,iy,iz] -= qz(DMgMg,CMg,ix,iy,iz) * Δzad_² +
                                qz(DMgFe,CFe,ix,iy,iz) * Δzad_² +
                                qz(DMgMn,CMn,ix,iy,iz) * Δzad_²
                @inbounds dtCFe[ix,iy,iz] -= qz(DFeMg,CMg,ix,iy,iz) * Δzad_² +
                                qz(DFeFe,CFe,ix,iy,iz) * Δzad_² +
                                qz(DFeMn,CMn,ix,iy,iz) * Δzad_²
                @inbounds dtCMn[ix,iy,iz] -= qz(DMnMg,CMg,ix,iy,iz) * Δzad_² +
                                qz(DMnFe,CFe,ix,iy,iz) * Δzad_² +
                                qz(DMnMn,CMn,ix,iy,iz) * Δzad_²
            end

        # if point is an inclusion or matrix
        else
            @inbounds dtCMg[ix,iy,iz] = 0.0
            @inbounds dtCFe[ix,iy,iz] = 0.0
            @inbounds dtCMn[ix,iy,iz] = 0.0
        end
        # if point is on grain boundary
        if isone(Grt_boundaries[ix,iy,iz])
            @inbounds dtCMg[ix,iy,iz] = 0.0
            @inbounds dtCFe[ix,iy,iz] = 0.0
            @inbounds dtCMn[ix,iy,iz] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        @inbounds dtCMg[ix,iy,iz] = 0.0
        @inbounds dtCFe[ix,iy,iz] = 0.0
        @inbounds dtCMn[ix,iy,iz] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian(du::Array_T,u::Array_T,p,t) where Array_T <: AbstractArray{<:Real, 4}

    (; D, D0, D_charact, Δxad_, Δyad_, Δzad_, diffcoef, D0_data, T, P, fugacity_O2, time_update_ad) = p.domain
    (; grt_position, grt_boundary) = p.domain.IC

    CMg = @view u[:,:,:,1]
    CFe = @view u[:,:,:,2]
    CMn = @view u[:,:,:,3]

    dtCMg = @view du[:,:,:,1]
    dtCFe = @view du[:,:,:,2]
    dtCMn = @view du[:,:,:,3]

    # find current T, P, and fugacity_O2
    # to do this, we need to find the lower bound of t in time_update_ad
    index = findfirst(x -> x >= t, time_update_ad)
    if index === nothing
        index = length(time_update_ad)
    end

    P_kbar = P[index] * 1u"kbar"
    T_K = (T[index]+273.15) * 1u"K"
    fO2 = (fugacity_O2[index])NoUnits

    D_charact_ = inv(D_charact)
    Δxad_² = Δxad_ * Δxad_
    Δyad_² = Δyad_ * Δyad_
    Δzad_² = Δzad_ * Δzad_

    # update diffusive parameters
    @parallel Diffusion_coef_3D_major!(D, CMg, CFe ,CMn, D0, D_charact_, grt_position, diffcoef, D0_data, T_K, P_kbar, fO2)

    # semi-discretization
    @parallel stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, grt_position, grt_boundary, Δxad_², Δyad_², Δzad_²)
end
