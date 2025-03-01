import Base.@propagate_inbounds

function Diffusion_coef_spherical_major!(D, CMg, CFe, CMn, D0, D_charact)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    D_charact_ = 1 / D_charact

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, I) = D0[1] * CMg[I] + D0[2] * CFe[I] + D0[3] * CMn[I] +
    D0[4] * (1 - CMg[I] - CFe[I] - CMn[I])

    @inbounds for I in eachindex(DMgMg)
        sum_D_ = 1 / (sum_D(CMg, CFe, CMn, D0, I))

        DMgMg[I] = (D0[1] - D0[1] * CMg[I] * sum_D_ * (D0[1] - D0[end])) * D_charact_
        DMgFe[I] = (      - D0[1] * CMg[I] * sum_D_ * (D0[2] - D0[end])) * D_charact_
        DMgMn[I] = (      - D0[1] * CMg[I] * sum_D_ * (D0[3] - D0[end])) * D_charact_

        DFeMg[I] = (      - D0[2] * CFe[I] * sum_D_ * (D0[1] - D0[end])) * D_charact_
        DFeFe[I] = (D0[2] - D0[2] * CFe[I] * sum_D_ * (D0[2] - D0[end])) * D_charact_
        DFeMn[I] = (      - D0[2] * CFe[I] * sum_D_ * (D0[3] - D0[end])) * D_charact_

        DMnMg[I] = (      - D0[3] * CMn[I] * sum_D_ * (D0[1] - D0[end])) * D_charact_
        DMnFe[I] = (      - D0[3] * CMn[I] * sum_D_ * (D0[2] - D0[end])) * D_charact_
        DMnMn[I] = (D0[3] - D0[3] * CMn[I] * sum_D_ * (D0[3] - D0[end])) * D_charact_
    end
end



function stencil_diffusion_spherical_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δrad, Δrad_, r_ad)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δrad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * Δrad_

    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, Δrad_)
        dtC[ix] = (qx(D1,C1,ix,Δrad_) - qx(D1,C1,ix-1,Δrad_)) * Δrad_ +
                  (qx(D2,C2,ix,Δrad_) - qx(D2,C2,ix-1,Δrad_)) * Δrad_ +
                  (qx(D3,C3,ix,Δrad_) - qx(D3,C3,ix-1,Δrad_)) * Δrad_ +
                  D1[ix] / r_ad[ix] * (C1[ix+1]-C1[ix-1]) * Δrad_ +
                  D2[ix] / r_ad[ix] * (C2[ix+1]-C2[ix-1]) * Δrad_ +
                  D3[ix] / r_ad[ix] * (C3[ix+1]-C3[ix-1]) * Δrad_
    end

    @inbounds for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, Δrad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, Δrad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, Δrad_)
        end
        # solve singularities (equivalent to homogeneous Neumann BC), see Versypt, A. N. F., & Braatz, R. D. (2014). Analysis of finite difference discretization schemes for diffusion in spheres with variable diffusivity. Computers & chemical engineering, 71, 241-252.
        if ix == 1
            dtCMg[ix] = (6 * DMgMg[ix] * (CMg[ix+1]-CMg[ix])) / (Δrad^2) +
                        (6 * DMgFe[ix] * (CFe[ix+1]-CFe[ix])) / (Δrad^2) +
                        (6 * DMgMn[ix] * (CMn[ix+1]-CMn[ix])) / (Δrad^2)
            dtCFe[ix] = (6 * DFeMg[ix] * (CMg[ix+1]-CMg[ix])) / (Δrad^2) +
                        (6 * DFeFe[ix] * (CFe[ix+1]-CFe[ix])) / (Δrad^2) +
                        (6 * DFeMn[ix] * (CMn[ix+1]-CMn[ix])) / (Δrad^2)
            dtCMn[ix] = (6 * DMnMg[ix] * (CMg[ix+1]-CMg[ix])) / (Δrad^2) +
                        (6 * DMnFe[ix] * (CFe[ix+1]-CFe[ix])) / (Δrad^2) +
                        (6 * DMnMn[ix] * (CMn[ix+1]-CMn[ix])) / (Δrad^2)
        end
    end
end

function semi_discretisation_diffusion_spherical(du,u,p,t)

    @unpack D, D0, D_charact, Δrad, Δrad_, r_ad = p.domain

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    Diffusion_coef_spherical_major!(D, CMg, CFe ,CMn, D0, D_charact)

    # semi-discretization
    stencil_diffusion_spherical_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δrad, Δrad_, r_ad)
end

