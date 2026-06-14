import Base.@propagate_inbounds

@parallel_indices (ix, iy, iz) function stencil_diffusion_3D_major!(
    dtCMg, dtCFe, dtCMn, CMg, CFe, CMn, D0, D_charact_,
    position_Grt, Grt_boundaries, diffcoef, D0_data, T_K, P_kbar, fO2,
    Δxad_², Δyad_², Δzad_²)

    # Compute the 9 interdiffusion coefficients from face-average composition.
    @propagate_inbounds @inline function interdiff(cMg, cFe, cMn, D0Mg, D0Fe, D0Mn, D0Ca)
        sum_D  = D0Mg * cMg + D0Fe * cFe + D0Mn * cMn + D0Ca * (1 - cMg - cFe - cMn)
        sum_D_ = inv(sum_D)
        sMg    = sum_D_ * (D0Mg - D0Ca)
        sFe    = sum_D_ * (D0Fe - D0Ca)
        sMn    = sum_D_ * (D0Mn - D0Ca)
        fMg    = D0Mg * cMg
        fFe    = D0Fe * cFe
        fMn    = D0Mn * cMn
        return (D0Mg - fMg * sMg) * D_charact_,
                   (-fMg * sFe)   * D_charact_,
                   (-fMg * sMn)   * D_charact_,
               (-fFe * sMg)       * D_charact_,
               (D0Fe - fFe * sFe) * D_charact_,
                   (-fFe * sMn)   * D_charact_,
               (-fMn * sMg)       * D_charact_,
                   (-fMn * sFe)   * D_charact_,
               (D0Mn - fMn * sMn) * D_charact_
    end

    # Compute flux at x-face between (ix,iy,iz) and (ix+1,iy,iz).
    @propagate_inbounds @inline function flux_x(ix, iy, iz)
        cMg = (CMg[ix,iy,iz] + CMg[ix+1,iy,iz]) / 2
        cFe = (CFe[ix,iy,iz] + CFe[ix+1,iy,iz]) / 2
        cMn = (CMn[ix,iy,iz] + CMn[ix+1,iy,iz]) / 2

        D0Mg_ = D0[1]; D0Fe_ = D0[2]; D0Mn_ = D0[3]; D0Ca_ = D0[4]
        if diffcoef == 2 || diffcoef == 3
            a0_Fe = eltype(CMg)(1.1525)
            a0_Mg = eltype(CMg)(1.1456)
            a0_Mn = eltype(CMg)(1.1614)
            a0_Ca = eltype(CMg)(1.1852)
            X = convert(eltype(CMg), cFe * a0_Fe + cMg * a0_Mg + cMn * a0_Mn + (1 - cMg - cFe - cMn) * a0_Ca)NoUnits
            D0Mg_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mg, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Fe_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Fe, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Mn_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mn, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Ca_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Ca, T=T_K, P=P_kbar, fO2=fO2, X=X))))
        end

        DMgMg, DMgFe, DMgMn,
        DFeMg, DFeFe, DFeMn,
        DMnMg, DMnFe, DMnMn = interdiff(cMg, cFe, cMn, D0Mg_, D0Fe_, D0Mn_, D0Ca_)
        ΔMg = CMg[ix+1,iy,iz] - CMg[ix,iy,iz]
        ΔFe = CFe[ix+1,iy,iz] - CFe[ix,iy,iz]
        ΔMn = CMn[ix+1,iy,iz] - CMn[ix,iy,iz]
        return DMgMg*ΔMg + DMgFe*ΔFe + DMgMn*ΔMn,
               DFeMg*ΔMg + DFeFe*ΔFe + DFeMn*ΔMn,
               DMnMg*ΔMg + DMnFe*ΔFe + DMnMn*ΔMn
    end

    # Compute flux at y-face between (ix,iy,iz) and (ix,iy+1,iz).
    @propagate_inbounds @inline function flux_y(ix, iy, iz)
        cMg = (CMg[ix,iy,iz] + CMg[ix,iy+1,iz]) / 2
        cFe = (CFe[ix,iy,iz] + CFe[ix,iy+1,iz]) / 2
        cMn = (CMn[ix,iy,iz] + CMn[ix,iy+1,iz]) / 2

        D0Mg_ = D0[1]; D0Fe_ = D0[2]; D0Mn_ = D0[3]; D0Ca_ = D0[4]
        if diffcoef == 2 || diffcoef == 3
            a0_Fe = eltype(CMg)(1.1525)
            a0_Mg = eltype(CMg)(1.1456)
            a0_Mn = eltype(CMg)(1.1614)
            a0_Ca = eltype(CMg)(1.1852)
            X = convert(eltype(CMg), cFe * a0_Fe + cMg * a0_Mg + cMn * a0_Mn + (1 - cMg - cFe - cMn) * a0_Ca)NoUnits
            D0Mg_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mg, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Fe_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Fe, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Mn_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mn, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Ca_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Ca, T=T_K, P=P_kbar, fO2=fO2, X=X))))
        end

        DMgMg, DMgFe, DMgMn,
        DFeMg, DFeFe, DFeMn,
        DMnMg, DMnFe, DMnMn = interdiff(cMg, cFe, cMn, D0Mg_, D0Fe_, D0Mn_, D0Ca_)
        ΔMg = CMg[ix,iy+1,iz] - CMg[ix,iy,iz]
        ΔFe = CFe[ix,iy+1,iz] - CFe[ix,iy,iz]
        ΔMn = CMn[ix,iy+1,iz] - CMn[ix,iy,iz]
        return DMgMg*ΔMg + DMgFe*ΔFe + DMgMn*ΔMn,
               DFeMg*ΔMg + DFeFe*ΔFe + DFeMn*ΔMn,
               DMnMg*ΔMg + DMnFe*ΔFe + DMnMn*ΔMn
    end

    # Compute flux at z-face between (ix,iy,iz) and (ix,iy,iz+1).
    @propagate_inbounds @inline function flux_z(ix, iy, iz)
        cMg = (CMg[ix,iy,iz] + CMg[ix,iy,iz+1]) / 2
        cFe = (CFe[ix,iy,iz] + CFe[ix,iy,iz+1]) / 2
        cMn = (CMn[ix,iy,iz] + CMn[ix,iy,iz+1]) / 2

        D0Mg_ = D0[1]; D0Fe_ = D0[2]; D0Mn_ = D0[3]; D0Ca_ = D0[4]
        if diffcoef == 2 || diffcoef == 3
            a0_Fe = eltype(CMg)(1.1525)
            a0_Mg = eltype(CMg)(1.1456)
            a0_Mn = eltype(CMg)(1.1614)
            a0_Ca = eltype(CMg)(1.1852)
            X = convert(eltype(CMg), cFe * a0_Fe + cMg * a0_Mg + cMn * a0_Mn + (1 - cMg - cFe - cMn) * a0_Ca)NoUnits
            D0Mg_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mg, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Fe_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Fe, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Mn_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Mn, T=T_K, P=P_kbar, fO2=fO2, X=X))))
            D0Ca_ = convert(eltype(CMg), ustrip(uconvert(u"µm^2/Myr", compute_D(D0_data.Grt_Ca, T=T_K, P=P_kbar, fO2=fO2, X=X))))
        end

        DMgMg, DMgFe, DMgMn,
        DFeMg, DFeFe, DFeMn,
        DMnMg, DMnFe, DMnMn = interdiff(cMg, cFe, cMn, D0Mg_, D0Fe_, D0Mn_, D0Ca_)
        ΔMg = CMg[ix,iy,iz+1] - CMg[ix,iy,iz]
        ΔFe = CFe[ix,iy,iz+1] - CFe[ix,iy,iz]
        ΔMn = CMn[ix,iy,iz+1] - CMn[ix,iy,iz]
        return DMgMg*ΔMg + DMgFe*ΔFe + DMgMn*ΔMn,
               DFeMg*ΔMg + DFeFe*ΔFe + DFeMn*ΔMn,
               DMnMg*ΔMg + DMnFe*ΔFe + DMnMn*ΔMn
    end

    if ix>1 && ix<size(dtCMg,1) && iy>1 && iy<size(dtCMg,2) && iz>1 && iz<size(dtCMg,3)
        if isone(position_Grt[ix,iy,iz]) && iszero(Grt_boundaries[ix,iy,iz])
            @inbounds begin
                fxr_Mg, fxr_Fe, fxr_Mn = flux_x(ix,   iy, iz)
                fxl_Mg, fxl_Fe, fxl_Mn = flux_x(ix-1, iy, iz)
                fyr_Mg, fyr_Fe, fyr_Mn = flux_y(ix, iy,   iz)
                fyl_Mg, fyl_Fe, fyl_Mn = flux_y(ix, iy-1, iz)
                fzr_Mg, fzr_Fe, fzr_Mn = flux_z(ix, iy, iz)
                fzl_Mg, fzl_Fe, fzl_Mn = flux_z(ix, iy, iz-1)

                dtCMg[ix,iy,iz] = (fxr_Mg - fxl_Mg) * Δxad_² + (fyr_Mg - fyl_Mg) * Δyad_² + (fzr_Mg - fzl_Mg) * Δzad_²
                dtCFe[ix,iy,iz] = (fxr_Fe - fxl_Fe) * Δxad_² + (fyr_Fe - fyl_Fe) * Δyad_² + (fzr_Fe - fzl_Fe) * Δzad_²
                dtCMn[ix,iy,iz] = (fxr_Mn - fxl_Mn) * Δxad_² + (fyr_Mn - fyl_Mn) * Δyad_² + (fzr_Mn - fzl_Mn) * Δzad_²

                # Zero-flux Neumann at garnet-matrix/inclusion interfaces.
                # Add back the face flux that was subtracted by the divergence.
                if iszero(position_Grt[ix-1,iy,iz])
                    dtCMg[ix,iy,iz] += fxl_Mg * Δxad_²
                    dtCFe[ix,iy,iz] += fxl_Fe * Δxad_²
                    dtCMn[ix,iy,iz] += fxl_Mn * Δxad_²
                end
                if iszero(position_Grt[ix+1,iy,iz])
                    dtCMg[ix,iy,iz] -= fxr_Mg * Δxad_²
                    dtCFe[ix,iy,iz] -= fxr_Fe * Δxad_²
                    dtCMn[ix,iy,iz] -= fxr_Mn * Δxad_²
                end
                if iszero(position_Grt[ix,iy-1,iz])
                    dtCMg[ix,iy,iz] += fyl_Mg * Δyad_²
                    dtCFe[ix,iy,iz] += fyl_Fe * Δyad_²
                    dtCMn[ix,iy,iz] += fyl_Mn * Δyad_²
                end
                if iszero(position_Grt[ix,iy+1,iz])
                    dtCMg[ix,iy,iz] -= fyr_Mg * Δyad_²
                    dtCFe[ix,iy,iz] -= fyr_Fe * Δyad_²
                    dtCMn[ix,iy,iz] -= fyr_Mn * Δyad_²
                end
                if iszero(position_Grt[ix,iy,iz-1])
                    dtCMg[ix,iy,iz] += fzl_Mg * Δzad_²
                    dtCFe[ix,iy,iz] += fzl_Fe * Δzad_²
                    dtCMn[ix,iy,iz] += fzl_Mn * Δzad_²
                end
                if iszero(position_Grt[ix,iy,iz+1])
                    dtCMg[ix,iy,iz] -= fzr_Mg * Δzad_²
                    dtCFe[ix,iy,iz] -= fzr_Fe * Δzad_²
                    dtCMn[ix,iy,iz] -= fzr_Mn * Δzad_²
                end
            end
        else
            @inbounds dtCMg[ix,iy,iz] = 0.0
            @inbounds dtCFe[ix,iy,iz] = 0.0
            @inbounds dtCMn[ix,iy,iz] = 0.0
        end
        if isone(Grt_boundaries[ix,iy,iz])
            @inbounds dtCMg[ix,iy,iz] = 0.0
            @inbounds dtCFe[ix,iy,iz] = 0.0
            @inbounds dtCMn[ix,iy,iz] = 0.0
        end
    else
        @inbounds dtCMg[ix,iy,iz] = 0.0
        @inbounds dtCFe[ix,iy,iz] = 0.0
        @inbounds dtCMn[ix,iy,iz] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian(du::Array_T, u::Array_T, p, t) where Array_T <: AbstractArray{<:Real, 4}

    (; D0, D_charact, Δxad_, Δyad_, Δzad_, diffcoef, D0_data, T, P, fugacity_O2, time_update_ad) = p.domain
    (; grt_position, grt_boundary) = p.domain.IC

    CMg = @view u[:,:,:,1]
    CFe = @view u[:,:,:,2]
    CMn = @view u[:,:,:,3]

    dtCMg = @view du[:,:,:,1]
    dtCFe = @view du[:,:,:,2]
    dtCMn = @view du[:,:,:,3]

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

    @parallel stencil_diffusion_3D_major!(
        dtCMg, dtCFe, dtCMn, CMg, CFe, CMn, D0, D_charact_,
        grt_position, grt_boundary, diffcoef, D0_data, T_K, P_kbar, fO2,
        Δxad_², Δyad_², Δzad_²)
end
