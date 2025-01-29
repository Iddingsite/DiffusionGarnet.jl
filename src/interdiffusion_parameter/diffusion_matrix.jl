
function D_ini!(D0,T,P, fugacity_O2=1e-25)

    Grt_Mg = Grt_Mg_Chakraborty1992
    Grt_Fe = Grt_Fe_Chakraborty1992
    Grt_Mn = Grt_Mn_Chakraborty1992

    Grt_Mg = SetChemicalDiffusion(Grt_Mg)
    Grt_Fe = SetChemicalDiffusion(Grt_Fe)
    Grt_Mn = SetChemicalDiffusion(Grt_Mn)

    T_K = (T+273.15) .* 1u"K"
    P_kbar = P * 1u"kbar"

    DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Mg, T = T_K, P = P_kbar)))
    DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Fe, T = T_K, P = P_kbar)))
    DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Mn, T = T_K, P = P_kbar)))
    DCa = 0.5 * DFe

    fugacity_ratio = fugacity_O2/1e-25  # current fO2 over fO2 buffered with graphite

    # after Chakraborty and Ganguly, 1991 (page 142, equation 2)
    DMg = exp(log(DMg) + 1/6 * log(fugacity_ratio))
    DFe = exp(log(DFe) + 1/6 * log(fugacity_ratio))
    DMn = exp(log(DMn) + 1/6 * log(fugacity_ratio))
    DCa = exp(log(DCa) + 1/6 * log(fugacity_ratio))

    D0 .= (DMg, DFe, DMn, DCa)   # in µm^2/Myr

end


function Diffusion_param!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, D0, CMg, CFe ,CMn)
    for I in CartesianIndices(DMgMg)
        DMgMg[I] = (D0[1] - D0[1] * CMg[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMgFe[I] = (      - D0[1] * CMg[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMgMn[I] = (      - D0[1] * CMg[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DFeMg[I] = (      - D0[2] * CFe[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DFeFe[I] = (D0[2] - D0[2] * CFe[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DFeMn[I] = (      - D0[2] * CFe[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DMnMg[I] = (      - D0[3] * CMn[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMnFe[I] = (      - D0[3] * CMn[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMnMn[I] = (D0[3] - D0[3] * CMn[I] / @sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact
    end
end

