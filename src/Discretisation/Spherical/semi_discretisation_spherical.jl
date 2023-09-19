import Base.@propagate_inbounds

function Diffusion_coef_spherical!(D, CMg, CFe, CMn, D0, D_charact)

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


function stencil_diffusion_spherical!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δrad, Δrad_, r_ad)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δrad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * Δrad_

    @inbounds for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            dtCMg[ix] = (qx(DMgMg,CMg,ix,Δrad_) - qx(DMgMg,CMg,ix-1,Δrad_)) * Δrad_ +
                        (qx(DMgFe,CFe,ix,Δrad_) - qx(DMgFe,CFe,ix-1,Δrad_)) * Δrad_ +
                        (qx(DMgMn,CMn,ix,Δrad_) - qx(DMgMn,CMn,ix-1,Δrad_)) * Δrad_ +
                        DMgMg[ix] / r_ad[ix] * (CMg[ix+1]-CMg[ix-1]) * Δrad_ +
                        DMgFe[ix] / r_ad[ix] * (CFe[ix+1]-CFe[ix-1]) * Δrad_ +
                        DMgMn[ix] / r_ad[ix] * (CMn[ix+1]-CMn[ix-1]) * Δrad_
            dtCFe[ix] = (qx(DFeMg,CMg,ix,Δrad_) - qx(DFeMg,CMg,ix-1,Δrad_)) * Δrad_ +
                        (qx(DFeFe,CFe,ix,Δrad_) - qx(DFeFe,CFe,ix-1,Δrad_)) * Δrad_ +
                        (qx(DFeMn,CMn,ix,Δrad_) - qx(DFeMn,CMn,ix-1,Δrad_)) * Δrad_ +
                        DFeMg[ix] / r_ad[ix] * (CMg[ix+1]-CMg[ix-1]) * Δrad_ +
                        DFeFe[ix] / r_ad[ix] * (CFe[ix+1]-CFe[ix-1]) * Δrad_ +
                        DFeMn[ix] / r_ad[ix] * (CMn[ix+1]-CMn[ix-1]) * Δrad_
            dtCMn[ix] = (qx(DMnMg,CMg,ix,Δrad_) - qx(DMnMg,CMg,ix-1,Δrad_)) * Δrad_ +
                        (qx(DMnFe,CFe,ix,Δrad_) - qx(DMnFe,CFe,ix-1,Δrad_)) * Δrad_ +
                        (qx(DMnMn,CMn,ix,Δrad_) - qx(DMnMn,CMn,ix-1,Δrad_)) * Δrad_ +
                        DMnMg[ix] / r_ad[ix] * (CMg[ix+1]-CMg[ix-1]) * Δrad_ +
                        DMnFe[ix] / r_ad[ix] * (CFe[ix+1]-CFe[ix-1]) * Δrad_ +
                        DMnMn[ix] / r_ad[ix] * (CMn[ix+1]-CMn[ix-1]) * Δrad_
        end

        # solve singularities
        if ix == 1
            dtCMg[ix] = 6 * DMgMg[ix] * (CMg[ix+1]-CMg[ix]) / (Δrad)^2 +
                        6 * DMgFe[ix] * (CFe[ix+1]-CFe[ix]) / (Δrad)^2 +
                        6 * DMgMn[ix] * (CMn[ix+1]-CMn[ix]) / (Δrad)^2
            dtCFe[ix] = 6 * DFeMg[ix] * (CMg[ix+1]-CMg[ix]) / (Δrad)^2 +
                        6 * DFeFe[ix] * (CFe[ix+1]-CFe[ix]) / (Δrad)^2 +
                        6 * DFeMn[ix] * (CMn[ix+1]-CMn[ix]) / (Δrad)^2
            dtCMn[ix] = 6 * DMnMg[ix] * (CMg[ix+1]-CMg[ix]) / (Δrad)^2 +
                        6 * DMnFe[ix] * (CFe[ix+1]-CFe[ix]) / (Δrad)^2 +
                        6 * DMnMn[ix] * (CMn[ix+1]-CMn[ix]) / (Δrad)^2
        end
    end
end

function semi_discretisation_diffusion_spherical(du,u,p,t)

    @unpack D, D0, D_charact, Δrad, Δrad_, r_ad = p

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    Diffusion_coef_1D!(D, CMg, CFe ,CMn, D0, D_charact)

    # semi-discretization
    stencil_diffusion_spherical!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δrad, Δrad_, r_ad)
end

