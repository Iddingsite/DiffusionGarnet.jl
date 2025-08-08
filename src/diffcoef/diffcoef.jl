

"""
    D_update!(D0, T_K, P_kbar, diffcoef, CMg, CFe, CMn, X, fugacity_O2=1e-25NoUnits)

Update the diffusion coefficients `D0` based on the temperature `T_K`, pressure `P_kbar`, and the chemical compositions `CMg`, `CFe`, and `CMn`. The optional parameter `fugacity_O2` is used to set the oxygen fugacity, defaulting to 1e-25 Pa (graphite buffer).
"""
function D_update!(D0, T_K, P_kbar, diffcoef, CMg, CFe, CMn, D0_data, fugacity_O2=1e-25NoUnits)  # by defaut 1e-25 Pa is graphite buffer


    X = 0NoUnits  # default value for CG92, no composition dependence

    if diffcoef == 2 || diffcoef == 3
        # end-member unit-cell dimensions
        a0_Fe = 1.1525
        a0_Mg = 1.1456
        a0_Mn = 1.1614
        a0_Ca = 1.1852

        X = (CFe * a0_Fe + CMg * a0_Mg + CMn * a0_Mn + (1 - (CMg + CFe + CMn)) * a0_Ca)NoUnits
    end

    DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
    DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
    DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))

    DCa = 0.0

    if diffcoef == 1
        DCa = 0.5 * DFe
    elseif diffcoef == 2 || diffcoef == 3
        DCa = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = (fugacity_O2)NoUnits, X = (X)NoUnits)))
    end

    D0 .= (float(DMg), float(DFe), float(DMn), float(DCa))   # in µm^2/Myr
end

@parallel_indices (ix, iy) function D_update_2D!(D0, T_K, P_kbar, diffcoef, CMg, CFe, CMn, D0_data, fugacity_O2, grt_position, grt_boundary)

    if ix>1 && ix<size(D0,2) && iy>1 && iy<size(D0,3)
        if grt_position[ix,iy] == 1.0 || grt_boundary[ix,iy] == 1.0

            # there is a composition dependence in the self-diffusion coefficients for C12 and CA15
            if diffcoef == 1

                X = 0NoUnits  # default value for CG92, no composition dependence

                D0_Mg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Fe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Mn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Ca = 0.5 * D0_Fe

                D0[1, ix, iy, iz] = D0_Mg
                D0[2, ix, iy, iz] = D0_Fe
                D0[3, ix, iy, iz] = D0_Mn
                D0[4, ix, iy, iz] = D0_Ca

            elseif diffcoef == 2 || diffcoef == 3

                a0_Fe = 1.1525
                a0_Mg = 1.1456
                a0_Mn = 1.1614
                a0_Ca = 1.1852

                X = CFe[ix,iy] * a0_Fe + CMg[ix,iy] * a0_Mg + CMn[ix,iy] * a0_Mn + (1 - (CMg[ix,iy] + CFe[ix,iy] + CMn[ix,iy])) * a0_Ca
                X = convert(Float64, X)NoUnits

                D0_Mg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Fe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Mn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Ca = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))

                D0[1, ix, iy, iz] = D0_Mg
                D0[2, ix, iy, iz] = D0_Fe
                D0[3, ix, iy, iz] = D0_Mn
                D0[4, ix, iy, iz] = D0_Ca
            end
        end
    end

    return
end

@parallel_indices (ix, iy, iz) function D_update_3D!(D0, T_K, P_kbar, diffcoef, CMg, CFe, CMn, D0_data, fugacity_O2, grt_position, grt_boundary)

    if ix>1 && ix<size(D0,2) && iy>1 && iy<size(D0,3) && iz>1 && iz<size(D0,4)
        if grt_position[ix,iy,iz] == 1.0 || grt_boundary[ix,iy,iz] == 1.0

            # there is a composition dependence in the self-diffusion coefficients for C12 and CA15
            if diffcoef == 1

                X = 0NoUnits  # default value for CG92, no composition dependence

                D0_Mg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Fe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Mn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Ca = 0.5 * D0_Fe

                D0[1, ix, iy, iz] = D0_Mg
                D0[2, ix, iy, iz] = D0_Fe
                D0[3, ix, iy, iz] = D0_Mn
                D0[4, ix, iy, iz] = D0_Ca

            elseif diffcoef == 2 || diffcoef == 3

                a0_Fe = 1.1525
                a0_Mg = 1.1456
                a0_Mn = 1.1614
                a0_Ca = 1.1852

                X = CFe[ix,iy] * a0_Fe + CMg[ix,iy] * a0_Mg + CMn[ix,iy] * a0_Mn + (1 - (CMg[ix,iy] + CFe[ix,iy] + CMn[ix,iy])) * a0_Ca
                X = convert(Float64, X)NoUnits
                D0_Mg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Fe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Mn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))
                D0_Ca = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = fugacity_O2, X = X)))

                D0[1, ix, iy, iz] = D0_Mg
                D0[2, ix, iy, iz] = D0_Fe
                D0[3, ix, iy, iz] = D0_Mn
                D0[4, ix, iy, iz] = D0_Ca
            end
        end
    end

    return
end