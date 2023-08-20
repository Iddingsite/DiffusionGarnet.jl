

function D_ini!(D0,T,P)

    R = 8.314462618  # constant of gas in J mol−1 K−1

    # Magnesium
    D0Mg = 1.1 * 1e-3 * 1e8  # pre-exponential constant in µm2 / yr
    Eₐ_Mg = 67997 * 4.1855 # activation energy at 1 bar in J / mol
    ΔV⁺Mg = 5.3  # activation volume in cm3 / mol

    # Iron
    D0Fe = 6.4 * 1e-4 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Fe = 65824 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Fe = 5.6  # activation volume in cm3 / mol

    # Manganese
    D0Mn = 5.1 * 1e-4 * 1e8  # pre-exponential constant in µm2 / yr
    Eₐ_Mn = 60569 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Mn = 6.0  # activation volume in cm3 / mol

    DMg = D0Mg * exp(- (Eₐ_Mg + (100 * (P-0.001) * ΔV⁺Mg)) / (R * T))  # in µm2 / s
    DFe = D0Fe * exp(- (Eₐ_Fe + (100 * (P-0.001) * ΔV⁺Fe)) / (R * T))  # in µm2 / s
    DMn = D0Mn * exp(- (Eₐ_Mn + (100 * (P-0.001) * ΔV⁺Mn)) / (R * T))  # in µm2 / s
    DCa = 0.5 * DFe

    D0 .= [DMg, DFe, DMn, DCa] .* (365.25 * 24 * 3600 * 1e6) # in Ma

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


@parallel_indices (iy, ix) function Diffusion_coef_2D!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, CMg, CFe ,CMn, DMg0, DFe0, DMn0, DCa0, position_Grt)

    Base.@propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, ix) = (D0[1] * CMg[ix] + D0[2] * CFe[ix] + D0[3] * CMn[ix] + D0[4] * (1 - CMg[ix] - CFe[ix] - CMn[ix]))

    if ix>1 && ix<size(DMgMg,2) && iy>1 && iy<size(DMgMg,1)
        if position_Grt[iy,ix] == 1.0
            DMgMg[iy,ix] = (DMg0 - DMg0 * CMg[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMg0 - DCa0)) / D_charact
            DMgFe[iy,ix] = (     - DMg0 * CMg[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DFe0 - DCa0)) / D_charact
            DMgMn[iy,ix] = (     - DMg0 * CMg[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMn0 - DCa0)) / D_charact
            DFeMg[iy,ix] = (     - DFe0 * CFe[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMg0 - DCa0)) / D_charact
            DFeFe[iy,ix] = (DFe0 - DFe0 * CFe[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DFe0 - DCa0)) / D_charact
            DFeMn[iy,ix] = (     - DFe0 * CFe[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMn0 - DCa0)) / D_charact
            DMnMg[iy,ix] = (     - DMn0 * CMn[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMg0 - DCa0)) / D_charact
            DMnFe[iy,ix] = (     - DMn0 * CMn[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DFe0 - DCa0)) / D_charact
            DMnMn[iy,ix] = (DMn0 - DMn0 * CMn[iy,ix] / @sum_D(CMg,CFe,CMn,DMg0,DFe0,DMn0,DCa0,ix,iy) *
                            (DMn0 - DCa0)) / D_charact
        else
            DMgMg[iy,ix] = 0
            DMgFe[iy,ix] = 0
            DMgMn[iy,ix] = 0
            DFeMg[iy,ix] = 0
            DFeFe[iy,ix] = 0
            DFeMn[iy,ix] = 0
            DMnMg[iy,ix] = 0
            DMnFe[iy,ix] = 0
            DMnMn[iy,ix] = 0
        end
    end
    return
end