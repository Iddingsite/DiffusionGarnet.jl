import Base.@propagate_inbounds

@parallel_indices (ix, iy, iz) function stencil_diffusion_3D_major!(dtC, C, D, position_Grt, Grt_boundaries, Δxad_, Δyad_, Δzad_)

    @propagate_inbounds @inline qx(D, C, ix, iy, iz, Δxad_) = D * (C[ix+1,iy,iz]-C[ix,iy,iz]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, iz, Δyad_) = D * (C[ix,iy+1,iz]-C[ix,iy,iz]) * Δyad_
    @propagate_inbounds @inline qz(D, C, ix, iy, iz, Δzad_) = D * (C[ix,iy,iz+1]-C[ix,iy,iz]) * Δzad_


    @propagate_inbounds @inline function update_dtC(dtC, D, C, ix, iy, iz, Δxad_, Δyad_, Δzad_)
        dtC[ix,iy,iz] = (qx(D,C,ix,iy,iz,Δxad_) - qx(D,C,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qy(D,C,ix,iy,iz,Δyad_) - qy(D,C,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qz(D,C,ix,iy,iz,Δzad_) - qz(D,C,ix,iy,iz-1,Δzad_)) * Δzad_
    end

    # iterate inside the arrays
    if ix>1 && ix<size(dtC,1) && iy>1 && iy<size(dtC,2) && iz>1 && iz<size(dtC,3)
        if position_Grt[ix,iy, iz] == 1.0 && Grt_boundaries[ix,iy, iz] == 0.0
            @inbounds update_dtC(dtC, D, C, ix, iy, iz, Δxad_, Δyad_, Δzad_)

            # first order Neumann if inclusions
            # south
            if position_Grt[ix-1,iy,iz] == 0.0
                @inbounds dtC[ix,iy,iz] -= - qx(D,C,ix-1,iy,iz,Δxad_) * Δxad_
            end
            # north
            if position_Grt[ix+1,iy,iz] == 0.0
                @inbounds dtC[ix,iy,iz] -= qx(D,C,ix,iy,iz,Δxad_) * Δxad_
            end
            # west
            if position_Grt[ix,iy-1,iz] == 0.0
                @inbounds dtC[ix,iy,iz] -= - qy(D,C,ix,iy-1,iz,Δyad_) * Δyad_
            end
            # east
            if position_Grt[ix,iy+1,iz] == 0.0
                @inbounds dtC[ix,iy,iz] -= qy(D,C,ix,iy,iz,Δyad_) * Δyad_
            end
            # bottom
            if position_Grt[ix,iy,iz-1] == 0.0
                @inbounds dtC[ix,iy,iz] -= - qz(D,C,ix,iy,iz-1,Δzad_) * Δzad_
            end
            # top
            if position_Grt[ix,iy,iz+1] == 0.0
                @inbounds dtC[ix,iy,iz] -= qz(D,C,ix,iy,iz,Δzad_) * Δzad_
            end

        # if point is an inclusion or matrix
        else
            @inbounds dtC[ix,iy,iz] = 0.0
        end
        # if point is on grain boundary
        if Grt_boundaries[ix,iy,iz] == 1.0
            @inbounds dtC[ix,iy,iz] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        @inbounds dtC[ix,iy,iz] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian_trace(du::Array_T,u::Array_T,p,t) where Array_T <: AbstractArray{<:Real, 3}

    @unpack D, D_charact, Δxad_, Δyad_, Δzad_ = p.domain
    @unpack grt_position, grt_boundary = p.domain.IC

    D_ad = D[1] / D_charact

    # semi-discretization
    @parallel stencil_diffusion_3D_major!(du, u, D_ad, grt_position, grt_boundary, Δxad_, Δyad_, Δzad_)
end
