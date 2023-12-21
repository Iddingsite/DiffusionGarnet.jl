
import Base.@propagate_inbounds

@parallel_indices (iy, ix) function Diffusion_coef_2D!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, CMg, CFe ,CMn, D0, D_charact, grt_position)

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, ix, iy) = D0[1] * CMg[iy, ix] + D0[2] * CFe[iy, ix] + D0[3] * CMn[iy, ix] +
        D0[4] * (1 - CMg[iy, ix] - CFe[iy, ix] - CMn[iy, ix])

    if ix>1 && ix<size(DMgMg,2) && iy>1 && iy<size(DMgMg,1)
        if grt_position[iy,ix] == 1.0
            DMgMg[iy,ix] = (D0[1] - D0[1] * CMg[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[1] - D0[end])) / D_charact
            DMgFe[iy,ix] = (      - D0[1] * CMg[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[2] - D0[end])) / D_charact
            DMgMn[iy,ix] = (      - D0[1] * CMg[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[3] - D0[end])) / D_charact
            DFeMg[iy,ix] = (      - D0[2] * CFe[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[1] - D0[end])) / D_charact
            DFeFe[iy,ix] = (D0[2] - D0[2] * CFe[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[2] - D0[end])) / D_charact
            DFeMn[iy,ix] = (      - D0[2] * CFe[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[3] - D0[end])) / D_charact
            DMnMg[iy,ix] = (      - D0[3] * CMn[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[1] - D0[end])) / D_charact
            DMnFe[iy,ix] = (      - D0[3] * CMn[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[2] - D0[end])) / D_charact
            DMnMn[iy,ix] = (D0[3] - D0[3] * CMn[iy,ix] / sum_D(CMg,CFe,CMn,D0,ix,iy) * (D0[3] - D0[end])) / D_charact
        else
            DMgMg[iy,ix] = 0.0
            DMgFe[iy,ix] = 0.0
            DMgMn[iy,ix] = 0.0
            DFeMg[iy,ix] = 0.0
            DFeFe[iy,ix] = 0.0
            DFeMn[iy,ix] = 0.0
            DMnMg[iy,ix] = 0.0
            DMnFe[iy,ix] = 0.0
            DMnMn[iy,ix] = 0.0
        end
    end

    return
end


@parallel_indices (iy, ix) function stencil_diffusion_2D!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, position_Grt, Grt_boundaries, Δxad_, Δyad_)

    @propagate_inbounds @inline av_D_x(D, ix, iy) = 0.5 * (D[iy,ix] + D[iy,ix+1])
    @propagate_inbounds @inline av_D_y(D, ix, iy) = 0.5 * (D[iy,ix] + D[iy+1,ix])
    @propagate_inbounds @inline qx(D, C, ix, iy, Δxad_) = av_D_x(D, ix, iy) * (C[iy,ix+1]-C[iy,ix]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, Δyad_) = av_D_y(D, ix, iy) * (C[iy+1,ix]-C[iy,ix]) * Δyad_


    # iterate inside the arrays
    if ix>1 && ix<size(dtCMg,2) && iy>1 && iy<size(dtCMg,1)
        if position_Grt[iy,ix] == 1.0
            dtCMg[iy,ix] = (qx(DMgMg,CMg,ix,iy,Δxad_) - qx(DMgMg,CMg,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DMgFe,CFe,ix,iy,Δxad_) - qx(DMgFe,CFe,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DMgMn,CMn,ix,iy,Δxad_) - qx(DMgMn,CMn,ix-1,iy,Δxad_)) * Δxad_ +
                           (qy(DMgMg,CMg,ix,iy,Δyad_) - qy(DMgMg,CMg,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DMgFe,CFe,ix,iy,Δyad_) - qy(DMgFe,CFe,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DMgMn,CMn,ix,iy,Δyad_) - qy(DMgMn,CMn,ix,iy-1,Δyad_)) * Δyad_
            dtCFe[iy,ix] = (qx(DFeMg,CMg,ix,iy,Δxad_) - qx(DFeMg,CMg,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DFeFe,CFe,ix,iy,Δxad_) - qx(DFeFe,CFe,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DFeMn,CMn,ix,iy,Δxad_) - qx(DFeMn,CMn,ix-1,iy,Δxad_)) * Δxad_ +
                           (qy(DFeMg,CMg,ix,iy,Δyad_) - qy(DFeMg,CMg,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DFeFe,CFe,ix,iy,Δyad_) - qy(DFeFe,CFe,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DFeMn,CMn,ix,iy,Δyad_) - qy(DFeMn,CMn,ix,iy-1,Δyad_)) * Δyad_
            dtCMn[iy,ix] = (qx(DMnMg,CMg,ix,iy,Δxad_) - qx(DMnMg,CMg,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DMnFe,CFe,ix,iy,Δxad_) - qx(DMnFe,CFe,ix-1,iy,Δxad_)) * Δxad_ +
                           (qx(DMnMn,CMn,ix,iy,Δxad_) - qx(DMnMn,CMn,ix-1,iy,Δxad_)) * Δxad_ +
                           (qy(DMnMg,CMg,ix,iy,Δyad_) - qy(DMnMg,CMg,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DMnFe,CFe,ix,iy,Δyad_) - qy(DMnFe,CFe,ix,iy-1,Δyad_)) * Δyad_ +
                           (qy(DMnMn,CMn,ix,iy,Δyad_) - qy(DMnMn,CMn,ix,iy-1,Δyad_)) * Δyad_
            # first order Neumann if inclusions
            # left
            if position_Grt[iy,ix-1] == 0.0
                dtCMg[iy,ix] -= - qx(DMgMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMgFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMgMn,CMn,ix-1,iy,Δxad_) * Δxad_
                dtCFe[iy,ix] -= - qx(DFeMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DFeFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DFeMn,CMn,ix-1,iy,Δxad_) * Δxad_
                dtCMn[iy,ix] -= - qx(DMnMg,CMg,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMnFe,CFe,ix-1,iy,Δxad_) * Δxad_ -
                                  qx(DMnMn,CMn,ix-1,iy,Δxad_) * Δxad_
            end
            # right
            if position_Grt[iy,ix+1] == 0.0
                dtCMg[iy,ix] -= qx(DMgMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DMgFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DMgMn,CMn,ix,iy,Δxad_) * Δxad_
                dtCFe[iy,ix] -= qx(DFeMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DFeFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DFeMn,CMn,ix,iy,Δxad_) * Δxad_
                dtCMn[iy,ix] -= qx(DMnMg,CMg,ix,iy,Δxad_) * Δxad_ +
                                qx(DMnFe,CFe,ix,iy,Δxad_) * Δxad_ +
                                qx(DMnMn,CMn,ix,iy,Δxad_) * Δxad_
            end
            # bottom
            if position_Grt[iy-1,ix] == 0.0
                dtCMg[iy,ix] -= - qy(DMgMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMgFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMgMn,CMn,ix,iy-1,Δyad_) * Δyad_
                dtCFe[iy,ix] -= - qy(DFeMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DFeFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DFeMn,CMn,ix,iy-1,Δyad_) * Δyad_
                dtCMn[iy,ix] -= - qy(DMnMg,CMg,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMnFe,CFe,ix,iy-1,Δyad_) * Δyad_ -
                                  qy(DMnMn,CMn,ix,iy-1,Δyad_) * Δyad_
            end
            # top
            if position_Grt[iy+1,ix] == 0.0
                dtCMg[iy,ix] -= qy(DMgMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DMgFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DMgMn,CMn,ix,iy,Δyad_) * Δyad_
                dtCFe[iy,ix] -= qy(DFeMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DFeFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DFeMn,CMn,ix,iy,Δyad_) * Δyad_
                dtCMn[iy,ix] -= qy(DMnMg,CMg,ix,iy,Δyad_) * Δyad_ +
                                qy(DMnFe,CFe,ix,iy,Δyad_) * Δyad_ +
                                qy(DMnMn,CMn,ix,iy,Δyad_) * Δyad_
            end
        # if point is an inclusion
        else
            dtCMg[iy,ix] = 0.0
            dtCFe[iy,ix] = 0.0
            dtCMn[iy,ix] = 0.0
        end
        # if point is on grain boundary
        if Grt_boundaries[iy,ix] == 1.0
            dtCMg[iy,ix] = 0.0
            dtCFe[iy,ix] = 0.0
            dtCMn[iy,ix] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        dtCMg[iy,ix] = 0.0
        dtCFe[iy,ix] = 0.0
        dtCMn[iy,ix] = 0.0
    end

    return
end


function semi_discretisation_diffusion_2D(du,u,p,t)

    @unpack D, D0, D_charact, Δxad_, Δyad_ = p.domain
    @unpack grt_position, grt_boundary = p.domain.IC
    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    CMg = @view u[:,:,1]
    CFe = @view u[:,:,2]
    CMn = @view u[:,:,3]

    dtCMg = @view du[:,:,1]
    dtCFe = @view du[:,:,2]
    dtCMn = @view du[:,:,3]

    # update diffusive parameters
    @parallel Diffusion_coef_2D!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 CMg, CFe ,CMn, D0, D_charact, grt_position)


    # semi-discretization
    @parallel stencil_diffusion_2D!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn,
                                 DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 grt_position, grt_boundary, Δxad_, Δyad_)
end

