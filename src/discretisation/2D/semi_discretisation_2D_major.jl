import Base.@propagate_inbounds

@parallel_indices (ix, iy) function Diffusion_coef_2D_major!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, CMg, CFe ,CMn, D0, D_charact, grt_position)

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, ix, iy) = D0[1] * CMg[ix, iy] + D0[2] * CFe[ix, iy] + D0[3] * CMn[ix, iy] +
        D0[4] * (1 - CMg[ix, iy] - CFe[ix, iy] - CMn[ix, iy])

    D_charact_ = 1 / D_charact

    if ix>1 && ix<size(DMgMg,1) && iy>1 && iy<size(DMgMg,2)
        if grt_position[ix,iy] == 1.0
            sum_D_ = 1 / (sum_D(CMg,CFe,CMn,D0,ix,iy))

            DMgMg[ix,iy] = (D0[1] - D0[1] * CMg[ix,iy] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DMgFe[ix,iy] = (      - D0[1] * CMg[ix,iy] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DMgMn[ix,iy] = (      - D0[1] * CMg[ix,iy] * sum_D_ * (D0[3] - D0[end])) * D_charact_
            DFeMg[ix,iy] = (      - D0[2] * CFe[ix,iy] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DFeFe[ix,iy] = (D0[2] - D0[2] * CFe[ix,iy] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DFeMn[ix,iy] = (      - D0[2] * CFe[ix,iy] * sum_D_ * (D0[3] - D0[end])) * D_charact_
            DMnMg[ix,iy] = (      - D0[3] * CMn[ix,iy] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DMnFe[ix,iy] = (      - D0[3] * CMn[ix,iy] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DMnMn[ix,iy] = (D0[3] - D0[3] * CMn[ix,iy] * sum_D_ * (D0[3] - D0[end])) * D_charact_
        end
    end

    return
end


@parallel_indices (ix, iy) function stencil_diffusion_2D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, position_Grt, Grt_boundaries, Δxad_, Δyad_)

    @propagate_inbounds @inline av_D_x(D, ix, iy)       = 0.5 * (D[ix,iy] + D[ix+1,iy])
    @propagate_inbounds @inline av_D_y(D, ix, iy)       = 0.5 * (D[ix,iy] + D[ix,iy+1])
    @propagate_inbounds @inline qx(D, C, ix, iy, Δxad_) = av_D_x(D, ix, iy) * (C[ix+1,iy]-C[ix,iy]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, Δyad_) = av_D_y(D, ix, iy) * (C[ix,iy+1]-C[ix,iy]) * Δyad_

    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, iy, Δxad_, Δyad_)
        dtC[ix,iy] = (qx(D1,C1,ix,iy,Δxad_) - qx(D1,C1,ix-1,iy,Δxad_)) * Δxad_ +
                     (qx(D2,C2,ix,iy,Δxad_) - qx(D2,C2,ix-1,iy,Δxad_)) * Δxad_ +
                     (qx(D3,C3,ix,iy,Δxad_) - qx(D3,C3,ix-1,iy,Δxad_)) * Δxad_ +
                     (qy(D1,C1,ix,iy,Δyad_) - qy(D1,C1,ix,iy-1,Δyad_)) * Δyad_ +
                     (qy(D2,C2,ix,iy,Δyad_) - qy(D2,C2,ix,iy-1,Δyad_)) * Δyad_ +
                     (qy(D3,C3,ix,iy,Δyad_) - qy(D3,C3,ix,iy-1,Δyad_)) * Δyad_
    end

    # iterate inside the arrays
    if ix>1 && ix<size(dtCMg,1) && iy>1 && iy<size(dtCMg,2)
        if position_Grt[ix,iy] == 1.0
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, iy, Δxad_, Δyad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, iy, Δxad_, Δyad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, iy, Δxad_, Δyad_)

            # first order Neumann if inclusions
            # bottom
            if position_Grt[ix-1,iy] == 0.0
                dtCMg[ix,iy] -= - qx(DMgMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMgFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMgMn,CMn,ix-1,iy,Δxad_) * Δxad_
                dtCFe[ix,iy] -= - qx(DFeMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DFeFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DFeMn,CMn,ix-1,iy,Δxad_) * Δxad_
                dtCMn[ix,iy] -= - qx(DMnMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMnFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMnMn,CMn,ix-1,iy,Δxad_) * Δxad_
            end
            # top
            if position_Grt[ix+1,iy] == 0.0
                dtCMg[ix,iy] -= qx(DMgMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DMgFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DMgMn,CMn,ix,iy,Δxad_) * Δxad_
                dtCFe[ix,iy] -= qx(DFeMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DFeFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DFeMn,CMn,ix,iy,Δxad_) * Δxad_
                dtCMn[ix,iy] -= qx(DMnMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DMnFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DMnMn,CMn,ix,iy,Δxad_) * Δxad_
            end
            # left
            if position_Grt[ix,iy-1] == 0.0
                dtCMg[ix,iy] -= - qy(DMgMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMgFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMgMn,CMn,ix,iy-1,Δyad_) * Δyad_
                dtCFe[ix,iy] -= - qy(DFeMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DFeFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DFeMn,CMn,ix,iy-1,Δyad_) * Δyad_
                dtCMn[ix,iy] -= - qy(DMnMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMnFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMnMn,CMn,ix,iy-1,Δyad_) * Δyad_
            end
            # right
            if position_Grt[ix,iy+1] == 0.0
                dtCMg[ix,iy] -= qy(DMgMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DMgFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DMgMn,CMn,ix,iy,Δyad_) * Δyad_
                dtCFe[ix,iy] -= qy(DFeMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DFeFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DFeMn,CMn,ix,iy,Δyad_) * Δyad_
                dtCMn[ix,iy] -= qy(DMnMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DMnFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DMnMn,CMn,ix,iy,Δyad_) * Δyad_
            end
        # if point is an inclusion or matrix
        else
            dtCMg[ix,iy] = 0.0
            dtCFe[ix,iy] = 0.0
            dtCMn[ix,iy] = 0.0
        end
        # if point is on grain boundary
        if Grt_boundaries[ix,iy] == 1.0
            dtCMg[ix,iy] = 0.0
            dtCFe[ix,iy] = 0.0
            dtCMn[ix,iy] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        dtCMg[ix,iy] = 0.0
        dtCFe[ix,iy] = 0.0
        dtCMn[ix,iy] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian(du::T,u::T,p,t) where T <: AbstractArray{<:Real, 3}

    @unpack D, D0, D_charact, Δxad_, Δyad_ = p.domain
    @unpack grt_position, grt_boundary     = p.domain.IC
    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    CMg = @view u[:,:,1]
    CFe = @view u[:,:,2]
    CMn = @view u[:,:,3]

    dtCMg = @view du[:,:,1]
    dtCFe = @view du[:,:,2]
    dtCMn = @view du[:,:,3]

    # update diffusive parameters
    @parallel Diffusion_coef_2D_major!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 CMg, CFe ,CMn, D0, D_charact, grt_position)


    # semi-discretization
    @parallel stencil_diffusion_2D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn,
                                 DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 grt_position, grt_boundary, Δxad_, Δyad_)
end

