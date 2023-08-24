import Base.@propagate_inbounds

function Diffusion_coef_1D!(D, CMg, CFe, CMn, D0, D_charact)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, I) = D0[1] * CMg[I] + D0[2] * CFe[I] + D0[3] * CMn[I] +
    D0[4] * (1 - CMg[I] - CFe[I] - CMn[I])

    for I in eachindex(DMgMg)
        DMgMg[I] = (D0[1] - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMgFe[I] = (      - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMgMn[I] = (      - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DFeMg[I] = (      - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DFeFe[I] = (D0[2] - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DFeMn[I] = (      - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DMnMg[I] = (      - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMnFe[I] = (      - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMnMn[I] = (D0[3] - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact
    end

end


function stencil_diffusion!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δxad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * Δxad_

    for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            dtCMg[ix] = (qx(DMgMg,CMg,ix,Δxad_) - qx(DMgMg,CMg,ix-1,Δxad_)) * Δxad_ +
                        (qx(DMgFe,CFe,ix,Δxad_) - qx(DMgFe,CFe,ix-1,Δxad_)) * Δxad_ +
                        (qx(DMgMn,CMn,ix,Δxad_) - qx(DMgMn,CMn,ix-1,Δxad_)) * Δxad_
            dtCFe[ix] = (qx(DFeMg,CMg,ix,Δxad_) - qx(DFeMg,CMg,ix-1,Δxad_)) * Δxad_ +
                        (qx(DFeFe,CFe,ix,Δxad_) - qx(DFeFe,CFe,ix-1,Δxad_)) * Δxad_ +
                        (qx(DFeMn,CMn,ix,Δxad_) - qx(DFeMn,CMn,ix-1,Δxad_)) * Δxad_
            dtCMn[ix] = (qx(DMnMg,CMg,ix,Δxad_) - qx(DMnMg,CMg,ix-1,Δxad_)) * Δxad_ +
                        (qx(DMnFe,CFe,ix,Δxad_) - qx(DMnFe,CFe,ix-1,Δxad_)) * Δxad_ +
                        (qx(DMnMn,CMn,ix,Δxad_) - qx(DMnMn,CMn,ix-1,Δxad_)) * Δxad_
        end
    end
end

function semi_dicretisation_diffusion_1D(du,u,p,t)

    @unpack D, D0, D_charact, Δxad_ = p

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    Diffusion_coef_1D!(D, CMg, CFe ,CMn, D0, D_charact)

    # semi-discretization
    stencil_diffusion!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_)
end

