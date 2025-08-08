

"""
    D_update!(D0, T, P, diffcoef, CMg, CFe, CMn, X, fugacity_O2=1e-25)

Update the diffusion coefficients `D0` based on the temperature `T`, pressure `P`, and the chemical compositions `CMg`, `CFe`, and `CMn`. The optional parameter `fugacity_O2` is used to set the oxygen fugacity, defaulting to 1e-25 Pa (graphite buffer).
"""
function D_update!(D0, T, P, diffcoef, CMg, CFe, CMn, D0_data, fugacity_O2=1e-25)  # by defaut 1e-25 Pa is graphite buffer

    T_K = (T+273.15) * 1u"K"
    P_kbar = P * 1u"kbar"

    X = 0NoUnits  # default value for CG92, no composition dependence

    if diffcoef == 2 || diffcoef == 3
        # end-member unit-cell dimensions
        a0_Fe = 1.1525
        a0_Mg = 1.1456
        a0_Mn = 1.1614
        a0_Ca = 1.1852

        X = (CFe * a0_Fe + CMg * a0_Mg + CMn * a0_Mn + (1 - (CMg + CFe + CMn)) * a0_Ca)NoUnits
    end

    DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = (fugacity_O2)NoUnits, X = X)))
    DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = (fugacity_O2)NoUnits, X = X)))
    DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = (fugacity_O2)NoUnits, X = X)))

    DCa = 0.0

    if diffcoef == 1
        DCa = 0.5 * DFe
    elseif diffcoef == 2 || diffcoef == 3
        DCa = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Ca, T = T_K, P = P_kbar, fO2 = (fugacity_O2)NoUnits, X = (X)NoUnits)))
    end

    D0 .= (float(DMg), float(DFe), float(DMn), float(DCa))   # in µm^2/Myr
end