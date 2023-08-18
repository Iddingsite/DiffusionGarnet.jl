


function Diffusion_param!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, D0, CMg, CFe ,CMn)
    for I in CartesianIndices(DMnMg)
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