import Base.@propagate_inbounds

@parallel_indices (ix, iy) function stencil_diffusion_2D_trace!(dtC, C, D, position_Grt, Grt_boundaries, Δxad_, Δyad_)

    @propagate_inbounds @inline qx(D, C, ix, iy, Δxad_) = D * (C[ix+1,iy]-C[ix,iy]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, Δyad_) = D * (C[ix,iy+1]-C[ix,iy]) * Δyad_

    @propagate_inbounds @inline function update_dtC!(dtC, D, C, ix, iy, Δxad_, Δyad_)
        dtC[ix,iy] = (qx(D,C,ix,iy,Δxad_) - qx(D,C,ix-1,iy,Δxad_)) * Δxad_ +
                     (qy(D,C,ix,iy,Δyad_) - qy(D,C,ix,iy-1,Δyad_)) * Δyad_
    end

    @inbounds begin
        # iterate inside the arrays
        if ix>1 && ix<size(dtC,1) && iy>1 && iy<size(dtC,2)
            if position_Grt[ix,iy] == 1.0
                update_dtC!(dtC, D, C, ix, iy, Δxad_, Δyad_)

                # first order Neumann if inclusions
                # bottom
                if position_Grt[ix-1,iy] == 0.0
                    dtC[ix,iy] -= - qx(D,C,ix-1,iy,Δxad_) * Δxad_
                end
                # top
                if position_Grt[ix+1,iy] == 0.0
                    dtC[ix,iy] -= qx(D,C,ix,iy,Δxad_) * Δxad_
                end
                # left
                if position_Grt[ix,iy-1] == 0.0
                    dtC[ix,iy] -= - qy(D,C,ix,iy-1,Δyad_) * Δyad_
                end
                # right
                if position_Grt[ix,iy+1] == 0.0
                    dtC[ix,iy] -= qy(D,C,ix,iy,Δyad_) * Δyad_
                end
            # if point is an inclusion or matrix
            else
                dtC[ix,iy] = 0.0
            end
            # if point is on grain boundary
            if Grt_boundaries[ix,iy] == 1.0
                dtC[ix,iy] = 0.0
            end
        # if point is on a model boundary (can define Neumann here)
        else
            dtC[ix,iy] = 0.0
        end
    end

    return
end

function semi_discretisation_diffusion_cartesian_trace(du::Array_T,u::Array_T,p,t) where Array_T <: AbstractArray{<:Real, 2}

    @unpack D, D_charact, Δxad_, Δyad_ = p.domain
    @unpack grt_position, grt_boundary = p.domain.IC

    D_ad = D / D_charact

    # semi-discretization
    @parallel stencil_diffusion_2D_trace!(du, u, D_ad, grt_position, grt_boundary, Δxad_, Δyad_)
end

