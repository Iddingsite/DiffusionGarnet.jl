
import Base.@propagate_inbounds

@parallel_indices (ix, iy, iz) function Diffusion_coef_3D_major!(D, CMg, CFe ,CMn, D0, D_charact_, grt_position, diffcoef, D0_data, T_K, P_kbar, fO2)

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0Mg, D0Fe, D0Mn, D0Ca, ix, iy, iz) = D0Mg * CMg[ix, iy, iz] + D0Fe * CFe[ix, iy, iz] + D0Mn * CMn[ix, iy, iz] +
        D0Ca * (1 - CMg[ix, iy, iz] - CFe[ix, iy, iz] - CMn[ix, iy, iz])

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    if ix>1 && ix<size(DMgMg,1) && iy>1 && iy<size(DMgMg,2) && iz>1 && iz<size(DMgMg,3)
        if grt_position[ix,iy,iz] == 1.0

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
                X = convert(Float64, X)NoUnits

                D0Mg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
                D0Fe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
                D0Mn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
                D0Ca = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = fO2, X = X)))
            end

            # factor that repeats 3 times for each diffusion coefficient
            sum_D_ = 1 / sum_D(CMg,CFe,CMn,D0Mg,D0Fe,D0Mn,D0Ca,ix,iy,iz)

            sum_DMg = sum_D_ * (D0Mg - D0Ca)
            sum_DFe = sum_D_ * (D0Fe - D0Ca)
            sum_DMn = sum_D_ * (D0Mn - D0Ca)

            # factor that repeats 3 times for each diffusion coefficient
            D0MgCMg = D0Mg * CMg[ix,iy,iz]
            D0FeCFe = D0Fe * CFe[ix,iy,iz]
            D0MnCMn = D0Mn * CMn[ix,iy,iz]

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

    return
end


@parallel_indices (ix, iy, iz) function stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, position_Grt, Grt_boundaries, Δxad_, Δyad_, Δzad_)

    @propagate_inbounds @inline av_D_x(D, ix, iy, iz) = 0.5 * (D[ix,iy,iz] + D[ix+1,iy,iz])
    @propagate_inbounds @inline av_D_y(D, ix, iy, iz) = 0.5 * (D[ix,iy,iz] + D[ix,iy+1,iz])
    @propagate_inbounds @inline av_D_z(D, ix, iy, iz) = 0.5 * (D[ix,iy,iz] + D[ix,iy,iz+1])
    @propagate_inbounds @inline qx(D, C, ix, iy, iz, Δxad_) = av_D_x(D, ix, iy, iz) * (C[ix+1,iy,iz]-C[ix,iy,iz]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, iz, Δyad_) = av_D_y(D, ix, iy, iz) * (C[ix,iy+1,iz]-C[ix,iy,iz]) * Δyad_
    @propagate_inbounds @inline qz(D, C, ix, iy, iz, Δzad_) = av_D_z(D, ix, iy, iz) * (C[ix,iy,iz+1]-C[ix,iy,iz]) * Δzad_


    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, iy, iz, Δxad_, Δyad_, Δzad_)
        dtC[ix,iy,iz] = (qx(D1,C1,ix,iy,iz,Δxad_) - qx(D1,C1,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qx(D2,C2,ix,iy,iz,Δxad_) - qx(D2,C2,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qx(D3,C3,ix,iy,iz,Δxad_) - qx(D3,C3,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qy(D1,C1,ix,iy,iz,Δyad_) - qy(D1,C1,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qy(D2,C2,ix,iy,iz,Δyad_) - qy(D2,C2,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qy(D3,C3,ix,iy,iz,Δyad_) - qy(D3,C3,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qz(D1,C1,ix,iy,iz,Δzad_) - qz(D1,C1,ix,iy,iz-1,Δzad_)) * Δzad_ +
                        (qz(D2,C2,ix,iy,iz,Δzad_) - qz(D2,C2,ix,iy,iz-1,Δzad_)) * Δzad_ +
                        (qz(D3,C3,ix,iy,iz,Δzad_) - qz(D3,C3,ix,iy,iz-1,Δzad_)) * Δzad_
    end

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    # iterate inside the arrays
    if ix>1 && ix<size(dtCMg,1) && iy>1 && iy<size(dtCMg,2) && iz>1 && iz<size(dtCMg,3)
        if position_Grt[ix,iy, iz] == 1.0 && Grt_boundaries[ix,iy, iz] == 0.0
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)

            # first order Neumann if inclusions
            # south
            if position_Grt[ix-1,iy,iz] == 0.0
                dtCMg[ix,iy,iz] -= - qx(DMgMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMgFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMgMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
                dtCFe[ix,iy,iz] -= - qx(DFeMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DFeFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DFeMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
                dtCMn[ix,iy,iz] -= - qx(DMnMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMnFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMnMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
            end
            # north
            if position_Grt[ix+1,iy,iz] == 0.0
                dtCMg[ix,iy,iz] -= qx(DMgMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMgFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMgMn,CMn,ix,iy,iz,Δxad_) * Δxad_
                dtCFe[ix,iy,iz] -= qx(DFeMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DFeFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DFeMn,CMn,ix,iy,iz,Δxad_) * Δxad_
                dtCMn[ix,iy,iz] -= qx(DMnMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMnFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMnMn,CMn,ix,iy,iz,Δxad_) * Δxad_
            end
            # west
            if position_Grt[ix,iy-1,iz] == 0.0
                dtCMg[ix,iy,iz] -= - qy(DMgMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMgFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMgMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
                dtCFe[ix,iy,iz] -= - qy(DFeMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DFeFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DFeMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
                dtCMn[ix,iy,iz] -= - qy(DMnMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMnFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMnMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
            end
            # east
            if position_Grt[ix,iy+1,iz] == 0.0
                dtCMg[ix,iy,iz] -= qy(DMgMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMgFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMgMn,CMn,ix,iy,iz,Δyad_) * Δyad_
                dtCFe[ix,iy,iz] -= qy(DFeMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DFeFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DFeMn,CMn,ix,iy,iz,Δyad_) * Δyad_
                dtCMn[ix,iy,iz] -= qy(DMnMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMnFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMnMn,CMn,ix,iy,iz,Δyad_) * Δyad_
            end
            # bottom
            if position_Grt[ix,iy,iz-1] == 0.0
                dtCMg[ix,iy,iz] -= - qz(DMgMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMgFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMgMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
                dtCFe[ix,iy,iz] -= - qz(DFeMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DFeFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DFeMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
                dtCMn[ix,iy,iz] -= - qz(DMnMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMnFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMnMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
            end
            # top
            if position_Grt[ix,iy,iz+1] == 0.0
                dtCMg[ix,iy,iz] -= qz(DMgMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMgFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMgMn,CMn,ix,iy,iz,Δzad_) * Δzad_
                dtCFe[ix,iy,iz] -= qz(DFeMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DFeFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DFeMn,CMn,ix,iy,iz,Δzad_) * Δzad_
                dtCMn[ix,iy,iz] -= qz(DMnMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMnFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMnMn,CMn,ix,iy,iz,Δzad_) * Δzad_
            end

        # if point is an inclusion or matrix
        else
            dtCMg[ix,iy,iz] = 0.0
            dtCFe[ix,iy,iz] = 0.0
            dtCMn[ix,iy,iz] = 0.0
        end
        # if point is on grain boundary
        if Grt_boundaries[ix,iy,iz] == 1.0
            dtCMg[ix,iy,iz] = 0.0
            dtCFe[ix,iy,iz] = 0.0
            dtCMn[ix,iy,iz] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        dtCMg[ix,iy,iz] = 0.0
        dtCFe[ix,iy,iz] = 0.0
        dtCMn[ix,iy,iz] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian(du::Array_T,u::Array_T,p,t) where Array_T <: AbstractArray{<:Real, 4}

    @unpack D, D0, D_charact, Δxad_, Δyad_, Δzad_, diffcoef, D0_data, T, P, fugacity_O2, time_update_ad = p.domain
    @unpack grt_position, grt_boundary = p.domain.IC

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

    D_charact_ = 1 / D_charact

    # update diffusive parameters
    @parallel Diffusion_coef_3D_major!(D, CMg, CFe ,CMn, D0, D_charact_, grt_position, diffcoef, D0_data, T_K, P_kbar, fO2)

    # semi-discretization
    @parallel stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, grt_position, grt_boundary, Δxad_, Δyad_, Δzad_)
end

