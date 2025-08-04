import Base.@propagate_inbounds

function stencil_diffusion_spherical_trace!(dtC, C, D, Δrad, Δrad_, r_ad)

    @propagate_inbounds @inline qx(D, C, ix, Δrad_) = D * (C[ix+1]-C[ix]) * Δrad_

    @propagate_inbounds @inline function update_dtC(dtC, D, C, ix, Δrad_)
        dtC[ix] = (qx(D,C,ix,Δrad_) - qx(D,C,ix-1,Δrad_)) * Δrad_ +
                  D[ix] / r_ad[ix] * (C[ix+1]-1[ix-1]) * Δrad_
    end

    @inbounds for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            update_dtC(dtC, D, C, ix, Δrad_)
        end
        # solve singularities (equivalent to homogeneous Neumann BC), see Versypt, A. N. F., & Braatz, R. D. (2014). Analysis of finite difference discretization schemes for diffusion in spheres with variable diffusivity. Computers & chemical engineering, 71, 241-252.
        if ix == 1
            dtCMg[ix] = (6 * D * (C[ix+1]-C[ix])) / (Δrad^2)
        end
    end
end

function semi_discretisation_diffusion_spherical(du,u,p,t)

    @unpack D, D_charact, Δrad, Δrad_, r_ad = p.domain

    D_ad = D / D_charact

    # semi-discretization
    stencil_diffusion_spherical_trace!(dt, C, D_ad, Δrad, Δrad_, r_ad)
end

