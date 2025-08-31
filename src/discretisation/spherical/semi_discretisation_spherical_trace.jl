import Base.@propagate_inbounds

function stencil_diffusion_spherical_trace!(dtC, C, D, Δr_ad_, r_ad)

    @propagate_inbounds @inline qx(D, C, ix, Δr_ad_) = D * (C[ix+1]-C[ix]) * Δr_ad_[ix]

    @propagate_inbounds @inline function update_dtC(dtC, D, C, ix, Δr_ad_)
        Δr_ad_centered = (Δr_ad_[ix] + Δr_ad_[ix-1]) / 2

        dtC[ix] = @muladd (qx(D,C,ix,Δr_ad_) - qx(D,C,ix-1,Δr_ad_)) * Δr_ad_centered +
                  2 * D / r_ad[ix] * (C[ix+1]-C[ix-1]) / (r_ad[ix+1] - r_ad[ix-1])
    end

    @inbounds for ix in eachindex(dtC)
        if ix > 1 && ix < size(dtC, 1)
            update_dtC(dtC, D, C, ix, Δr_ad_)
        end
        # solve singularities (equivalent to homogeneous Neumann BC), see Versypt, A. N. F., & Braatz, R. D. (2014). Analysis of finite difference discretization schemes for diffusion in spheres with variable diffusivity. Computers & chemical engineering, 71, 241-252.
        if ix == 1
            if r_ad[1] == 0.0
                dtC[ix] = @muladd (6 * D * (C[ix+1]-C[ix])) * (Δr_ad_[1]^2)
            end
        end
    end
end

function semi_discretisation_diffusion_spherical_trace(du::T,u::T,p,t) where T <: AbstractArray{<:Real, 1}

    @unpack D, D_charact, Δr_ad, Δr_ad_, r_ad = p.domain

    D_ad = D / D_charact

    # semi-discretization
    stencil_diffusion_spherical_trace!(du, u, D_ad, Δr_ad_, r_ad)
end

