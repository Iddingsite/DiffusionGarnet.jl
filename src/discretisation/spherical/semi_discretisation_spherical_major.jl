import Base.@propagate_inbounds

function stencil_diffusion_spherical_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δr_ad_, r_ad)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δr_ad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * Δr_ad_[ix]

    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, Δr_ad_)
        Δr_ad_centered = (Δr_ad_[ix] + Δr_ad_[ix-1]) / 2

        dtC[ix] = (qx(D1,C1,ix,Δr_ad_) - qx(D1,C1,ix-1,Δr_ad_)) * Δr_ad_centered +
                  (qx(D2,C2,ix,Δr_ad_) - qx(D2,C2,ix-1,Δr_ad_)) * Δr_ad_centered +
                  (qx(D3,C3,ix,Δr_ad_) - qx(D3,C3,ix-1,Δr_ad_)) * Δr_ad_centered +
                  2 * D1[ix] / r_ad[ix] * (C1[ix+1]-C1[ix-1]) / (r_ad[ix+1] - r_ad[ix-1]) +
                  2 * D2[ix] / r_ad[ix] * (C2[ix+1]-C2[ix-1]) / (r_ad[ix+1] - r_ad[ix-1]) +
                  2 * D3[ix] / r_ad[ix] * (C3[ix+1]-C3[ix-1]) / (r_ad[ix+1] - r_ad[ix-1])
    end

    for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, Δr_ad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, Δr_ad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, Δr_ad_)
        end

        if ix == 1

            # solve singularities (equivalent to homogeneous Neumann BC), see Versypt, A. N. F., & Braatz, R. D. (2014). Analysis of finite difference discretization schemes for diffusion in spheres with variable diffusivity. Computers & chemical engineering, 71, 241-252.
            # only if points start at 0
            if r_ad[1] == 0.0
                dtCMg[ix] = (6 * DMgMg[ix] * (CMg[ix+1]-CMg[ix])) * (Δr_ad_[1]^2) +
                            (6 * DMgFe[ix] * (CFe[ix+1]-CFe[ix])) * (Δr_ad_[1]^2) +
                            (6 * DMgMn[ix] * (CMn[ix+1]-CMn[ix])) * (Δr_ad_[1]^2)
                dtCFe[ix] = (6 * DFeMg[ix] * (CMg[ix+1]-CMg[ix])) * (Δr_ad_[1]^2) +
                            (6 * DFeFe[ix] * (CFe[ix+1]-CFe[ix])) * (Δr_ad_[1]^2) +
                            (6 * DFeMn[ix] * (CMn[ix+1]-CMn[ix])) * (Δr_ad_[1]^2)
                dtCMn[ix] = (6 * DMnMg[ix] * (CMg[ix+1]-CMg[ix])) * (Δr_ad_[1]^2) +
                            (6 * DMnFe[ix] * (CFe[ix+1]-CFe[ix])) * (Δr_ad_[1]^2) +
                            (6 * DMnMn[ix] * (CMn[ix+1]-CMn[ix])) * (Δr_ad_[1]^2)
            end
        end
    end
end

function semi_discretisation_diffusion_spherical(du,u,p,t)

    @unpack D, D0, D_charact, Δr_ad_, r_ad = p.domain

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters (we can use the same function as in 1D)
    Diffusion_coef_1D_major!(D, CMg, CFe, CMn, D0, D_charact, p.domain, t)

    # semi-discretization
    stencil_diffusion_spherical_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δr_ad_, r_ad)
end

