
import Base.@propagate_inbounds

@parallel_indices (iy, ix, iz) function Diffusion_coef_3D_major!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, CMg, CFe ,CMn, D0, D_charact, grt_position)

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, ix, iy, iz) = D0[1] * CMg[iy, ix, iz] + D0[2] * CFe[iy, ix, iz] + D0[3] * CMn[iy, ix, iz] +
        D0[4] * (1 - CMg[iy, ix, iz] - CFe[iy, ix, iz] - CMn[iy, ix, iz])

    D_charact_ = 1 / D_charact

    if ix>1 && ix<size(DMgMg,2) && iy>1 && iy<size(DMgMg,1) && iz>1 && iz<size(DMgMg,3)
        if grt_position[iy,ix,iz] == 1.0
            sum_D_ = 1 / sum_D(CMg,CFe,CMn,D0,ix,iy, iz)
            DMgMg[iy,ix,iz] = (D0[1] - D0[1] * CMg[iy,ix,iz] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DMgFe[iy,ix,iz] = (      - D0[1] * CMg[iy,ix,iz] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DMgMn[iy,ix,iz] = (      - D0[1] * CMg[iy,ix,iz] * sum_D_ * (D0[3] - D0[end])) * D_charact_
            DFeMg[iy,ix,iz] = (      - D0[2] * CFe[iy,ix,iz] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DFeFe[iy,ix,iz] = (D0[2] - D0[2] * CFe[iy,ix,iz] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DFeMn[iy,ix,iz] = (      - D0[2] * CFe[iy,ix,iz] * sum_D_ * (D0[3] - D0[end])) * D_charact_
            DMnMg[iy,ix,iz] = (      - D0[3] * CMn[iy,ix,iz] * sum_D_ * (D0[1] - D0[end])) * D_charact_
            DMnFe[iy,ix,iz] = (      - D0[3] * CMn[iy,ix,iz] * sum_D_ * (D0[2] - D0[end])) * D_charact_
            DMnMn[iy,ix,iz] = (D0[3] - D0[3] * CMn[iy,ix,iz] * sum_D_ * (D0[3] - D0[end])) * D_charact_
        end
    end

    return
end


@parallel_indices (iy, ix, iz) function stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn, position_Grt, Grt_boundaries, Δxad_, Δyad_, Δzad_)

    @propagate_inbounds @inline av_D_x(D, ix, iy, iz) = 0.5 * (D[iy,ix,iz] + D[iy,ix+1,iz])
    @propagate_inbounds @inline av_D_y(D, ix, iy, iz) = 0.5 * (D[iy,ix,iz] + D[iy+1,ix,iz])
    @propagate_inbounds @inline av_D_z(D, ix, iy, iz) = 0.5 * (D[iy,ix,iz] + D[iy,ix,iz+1])
    @propagate_inbounds @inline qx(D, C, ix, iy, iz, Δxad_) = av_D_x(D, ix, iy, iz) * (C[iy,ix+1,iz]-C[iy,ix,iz]) * Δxad_
    @propagate_inbounds @inline qy(D, C, ix, iy, iz, Δyad_) = av_D_y(D, ix, iy, iz) * (C[iy+1,ix,iz]-C[iy,ix,iz]) * Δyad_
    @propagate_inbounds @inline qz(D, C, ix, iy, iz, Δzad_) = av_D_z(D, ix, iy, iz) * (C[iy,ix,iz+1]-C[iy,ix,iz]) * Δzad_


    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, iy, iz, Δxad_, Δyad_, Δzad_)
        dtC[iy,ix,iz] = (qx(D1,C1,ix,iy,iz,Δxad_) - qx(D1,C1,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qx(D2,C2,ix,iy,iz,Δxad_) - qx(D2,C2,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qx(D3,C3,ix,iy,iz,Δxad_) - qx(D3,C3,ix-1,iy,iz,Δxad_)) * Δxad_ +
                        (qy(D1,C1,ix,iy,iz,Δyad_) - qy(D1,C1,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qy(D2,C2,ix,iy,iz,Δyad_) - qy(D2,C2,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qy(D3,C3,ix,iy,iz,Δyad_) - qy(D3,C3,ix,iy-1,iz,Δyad_)) * Δyad_ +
                        (qz(D1,C1,ix,iy,iz,Δzad_) - qz(D1,C1,ix,iy,iz-1,Δzad_)) * Δzad_ +
                        (qz(D2,C2,ix,iy,iz,Δzad_) - qz(D2,C2,ix,iy,iz-1,Δzad_)) * Δzad_ +
                        (qz(D3,C3,ix,iy,iz,Δzad_) - qz(D3,C3,ix,iy,iz-1,Δzad_)) * Δzad_
    end


    # iterate inside the arrays
    if ix>1 && ix<size(dtCMg,2) && iy>1 && iy<size(dtCMg,1) && iz>1 && iz<size(dtCMg,3)
        if position_Grt[iy,ix, iz] == 1.0 && Grt_boundaries[iy,ix, iz] == 0.0
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, iy, iz, Δxad_, Δyad_, Δzad_)

            # first order Neumann if inclusions
            # west
            if position_Grt[iy,ix-1,iz] == 0.0
                dtCMg[iy,ix,iz] -= - qx(DMgMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMgFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMgMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
                dtCFe[iy,ix,iz] -= - qx(DFeMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DFeFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DFeMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
                dtCMn[iy,ix,iz] -= - qx(DMnMg,CMg,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMnFe,CFe,ix-1,iy,iz,Δxad_) * Δxad_ -
                                  qx(DMnMn,CMn,ix-1,iy,iz,Δxad_) * Δxad_
            end
            # east
            if position_Grt[iy,ix+1,iz] == 0.0
                dtCMg[iy,ix,iz] -= qx(DMgMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMgFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMgMn,CMn,ix,iy,iz,Δxad_) * Δxad_
                dtCFe[iy,ix,iz] -= qx(DFeMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DFeFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DFeMn,CMn,ix,iy,iz,Δxad_) * Δxad_
                dtCMn[iy,ix,iz] -= qx(DMnMg,CMg,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMnFe,CFe,ix,iy,iz,Δxad_) * Δxad_ +
                                qx(DMnMn,CMn,ix,iy,iz,Δxad_) * Δxad_
            end
            # south
            if position_Grt[iy-1,ix,iz] == 0.0
                dtCMg[iy,ix,iz] -= - qy(DMgMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMgFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMgMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
                dtCFe[iy,ix,iz] -= - qy(DFeMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DFeFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DFeMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
                dtCMn[iy,ix,iz] -= - qy(DMnMg,CMg,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMnFe,CFe,ix,iy-1,iz,Δyad_) * Δyad_ -
                                  qy(DMnMn,CMn,ix,iy-1,iz,Δyad_) * Δyad_
            end
            # north
            if position_Grt[iy+1,ix,iz] == 0.0
                dtCMg[iy,ix,iz] -= qy(DMgMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMgFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMgMn,CMn,ix,iy,iz,Δyad_) * Δyad_
                dtCFe[iy,ix,iz] -= qy(DFeMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DFeFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DFeMn,CMn,ix,iy,iz,Δyad_) * Δyad_
                dtCMn[iy,ix,iz] -= qy(DMnMg,CMg,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMnFe,CFe,ix,iy,iz,Δyad_) * Δyad_ +
                                qy(DMnMn,CMn,ix,iy,iz,Δyad_) * Δyad_
            end
            # bottom
            if position_Grt[iy,ix,iz-1] == 0.0
                dtCMg[iy,ix,iz] -= - qz(DMgMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMgFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMgMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
                dtCFe[iy,ix,iz] -= - qz(DFeMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DFeFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DFeMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
                dtCMn[iy,ix,iz] -= - qz(DMnMg,CMg,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMnFe,CFe,ix,iy,iz-1,Δzad_) * Δzad_ -
                                  qz(DMnMn,CMn,ix,iy,iz-1,Δzad_) * Δzad_
            end
            # top
            if position_Grt[iy,ix,iz+1] == 0.0
                dtCMg[iy,ix,iz] -= qz(DMgMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMgFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMgMn,CMn,ix,iy,iz,Δzad_) * Δzad_
                dtCFe[iy,ix,iz] -= qz(DFeMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DFeFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DFeMn,CMn,ix,iy,iz,Δzad_) * Δzad_
                dtCMn[iy,ix,iz] -= qz(DMnMg,CMg,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMnFe,CFe,ix,iy,iz,Δzad_) * Δzad_ +
                                qz(DMnMn,CMn,ix,iy,iz,Δzad_) * Δzad_
            end

        # if point is an inclusion or matrix
        else
            dtCMg[iy,ix,iz] = 0.0
            dtCFe[iy,ix,iz] = 0.0
            dtCMn[iy,ix,iz] = 0.0
        end
        # if point is on grain boundary
        if Grt_boundaries[iy,ix,iz] == 1.0
            dtCMg[iy,ix,iz] = 0.0
            dtCFe[iy,ix,iz] = 0.0
            dtCMn[iy,ix,iz] = 0.0
        end
    # if point is on a model boundary (can define Neumann here)
    else
        dtCMg[iy,ix,iz] = 0.0
        dtCFe[iy,ix,iz] = 0.0
        dtCMn[iy,ix,iz] = 0.0
    end

    return
end


function semi_discretisation_diffusion_cartesian(du::T,u::T,p,t) where T <: AbstractArray{<:Real, 4}

    @unpack D, D0, D_charact, Δxad_, Δyad_, Δzad_ = p.domain
    @unpack grt_position, grt_boundary = p.domain.IC
    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    CMg = @view u[:,:,:,1]
    CFe = @view u[:,:,:,2]
    CMn = @view u[:,:,:,3]

    dtCMg = @view du[:,:,:,1]
    dtCFe = @view du[:,:,:,2]
    dtCMn = @view du[:,:,:,3]

    # update diffusive parameters
    @parallel Diffusion_coef_3D_major!(DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 CMg, CFe ,CMn, D0, D_charact, grt_position)


    # semi-discretization
    @parallel stencil_diffusion_3D_major!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn,
                                 DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn,
                                 grt_position, grt_boundary, Δxad_, Δyad_, Δzad_)
end

