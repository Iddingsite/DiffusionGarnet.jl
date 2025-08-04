import Base.@propagate_inbounds

function stencil_diffusion_1D_trace!(dtC, C, D, Δxad_, bc_neumann)

    @propagate_inbounds @inline qx(D, C, ix, Δxad_) = D * (C[ix+1]-C[ix]) * Δxad_
    @propagate_inbounds @inline neumann_left(D, C, ix, Δxad_) = - D * (C[ix]-C[ix+1]) * Δxad_
    @propagate_inbounds @inline neumann_right(D, C, ix, Δxad_) = D * (C[ix-1]-C[ix]) * Δxad_

    @propagate_inbounds @inline function update_dtC(dtC, D, C, ix, Δxad_)
        dtC[ix] = (qx(D,C,ix,Δxad_) - qx(D,C,ix-1,Δxad_)) * Δxad_
    end

    @inbounds for ix in eachindex(dtC)
        if ix > 1 && ix < size(dtC, 1)
            update_dtC(dtC, D, C, ix, Δxad_)
        end
    end

    # neumann boundary conditions
    if bc_neumann[1] == true
        ix = 1
        dtC[1] = qx(D,C,ix,Δxad_) * Δxad
        dtC[1] += neumann_left(D,C,ix,Δxad_) * Δxad_
    end

    if bc_neumann[2] == true
        ix = last_index(dtC)
        dtC[end] = - qx(D,C,ix,Δxad_) * Δxad_
        dtC[end] += neumann_right(D,C,ix,Δxad_) * Δxad_
    end
end


function semi_discretisation_diffusion_cartesian_trace(du::T,u::T,p,t) where T <: AbstractArray{<:Real, 1}

    @unpack D, Δxad_, bc_neumann, D_charact = p.domain

    D_ad = D / D_charact

    # semi-discretization
    stencil_diffusion_1D_trace!(du, u, D_ad, Δxad_, bc_neumann)
end

