import Base.@propagate_inbounds

function Diffusion_coef_1D!(D, CMg, CFe, CMn, D0, D_charact)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline sum_D(CMg, CFe, CMn, D0, I) = D0[1] * CMg[I] + D0[2] * CFe[I] + D0[3] * CMn[I] +
    D0[4] * (1 - CMg[I] - CFe[I] - CMn[I])

    for I in eachindex(DMgMg)
        DMgMg[I] = (D0[1] - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMgFe[I] = (      - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMgMn[I] = (      - D0[1] * CMg[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DFeMg[I] = (      - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DFeFe[I] = (D0[2] - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DFeMn[I] = (      - D0[2] * CFe[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact

        DMnMg[I] = (      - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[1] - D0[end])) / D_charact
        DMnFe[I] = (      - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[2] - D0[end])) / D_charact
        DMnMn[I] = (D0[3] - D0[3] * CMn[I] / sum_D(CMg, CFe, CMn, D0, I) * (D0[3] - D0[end])) / D_charact
    end

end


function stencil_diffusion_1D!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_, bc_neumann)

    DMgMg, DMgFe, DMgMn, DFeMg, DFeFe, DFeMn, DMnMg, DMnFe, DMnMn = D

    @propagate_inbounds @inline qx(D, C, ix, Δxad_) = 0.5 * (D[ix] + D[ix+1]) * (C[ix+1]-C[ix]) * Δxad_
    @propagate_inbounds @inline neumann_left(D, C, ix, Δxad_) = - D[ix] * (C[ix]-C[ix+1]) * Δxad_
    @propagate_inbounds @inline neumann_right(D, C, ix, Δxad_) = D[ix] * (C[ix-1]-C[ix]) * Δxad_

    @propagate_inbounds @inline function update_dtC(dtC, D1, D2, D3, C1, C2, C3, ix, Δxad_)
        dtC[ix] = (qx(D1,C1,ix,Δxad_) - qx(D1,C1,ix-1,Δxad_)) * Δxad_ +
                  (qx(D2,C2,ix,Δxad_) - qx(D2,C2,ix-1,Δxad_)) * Δxad_ +
                  (qx(D3,C3,ix,Δxad_) - qx(D3,C3,ix-1,Δxad_)) * Δxad_
    end

    @inbounds for ix in eachindex(dtCMg)
        if ix > 1 && ix < size(dtCMg, 1)
            update_dtC(dtCMg, DMgMg, DMgFe, DMgMn, CMg, CFe, CMn, ix, Δxad_)
            update_dtC(dtCFe, DFeMg, DFeFe, DFeMn, CMg, CFe, CMn, ix, Δxad_)
            update_dtC(dtCMn, DMnMg, DMnFe, DMnMn, CMg, CFe, CMn, ix, Δxad_)
        end
    end

    # neumann boundary conditions
    if bc_neumann[1] == true
        ix = 1
        dtCMg[1] = qx(DMgMg,CMg,ix,Δxad_) * Δxad_ +
                   qx(DMgFe,CFe,ix,Δxad_) * Δxad_ +
                   qx(DMgMn,CMn,ix,Δxad_) * Δxad_
        dtCFe[1] = qx(DFeMg,CMg,ix,Δxad_) * Δxad_ +
                   qx(DFeFe,CFe,ix,Δxad_) * Δxad_ +
                   qx(DFeMn,CMn,ix,Δxad_) * Δxad_
        dtCMn[1] = qx(DMnMg,CMg,ix,Δxad_) * Δxad_ +
                   qx(DMnFe,CFe,ix,Δxad_) * Δxad_ +
                   qx(DMnMn,CMn,ix,Δxad_) * Δxad_
        dtCMg[1] += neumann_left(DMgMg,CMg,ix,Δxad_) * Δxad_ +
                    neumann_left(DMgFe,CFe,ix,Δxad_) * Δxad_ +
                    neumann_left(DMgMn,CMn,ix,Δxad_) * Δxad_
        dtCFe[1] += neumann_left(DFeMg,CMg,ix,Δxad_) * Δxad_ +
                    neumann_left(DFeFe,CFe,ix,Δxad_) * Δxad_ +
                    neumann_left(DFeMn,CMn,ix,Δxad_) * Δxad_
        dtCMn[1] += neumann_left(DMnMg,CMg,ix,Δxad_) * Δxad_ +
                    neumann_left(DMnFe,CFe,ix,Δxad_) * Δxad_ +
                    neumann_left(DMnMn,CMn,ix,Δxad_) * Δxad_
    end

    if bc_neumann[2] == true
        ix = last_index(dtCMg)
        dtCMg[end] = - qx(DMgMg,CMg,ix,Δxad_) * Δxad_ -
                       qx(DMgFe,CFe,ix,Δxad_) * Δxad_ -
                       qx(DMgMn,CMn,ix,Δxad_) * Δxad_
        dtCFe[end] = - qx(DFeMg,CMg,ix,Δxad_) * Δxad_ -
                       qx(DFeFe,CFe,ix,Δxad_) * Δxad_ -
                       qx(DFeMn,CMn,ix,Δxad_) * Δxad_
        dtCMn[end] = - qx(DMnMg,CMg,ix,Δxad_) * Δxad_ -
                       qx(DMnFe,CFe,ix,Δxad_) * Δxad_ -
                       qx(DMnMn,CMn,ix,Δxad_) * Δxad_
        dtCMg[end] += neumann_right(DMgMg,CMg,ix,Δxad_) * Δxad_ +
                      neumann_right(DMgFe,CFe,ix,Δxad_) * Δxad_ +
                      neumann_right(DMgMn,CMn,ix,Δxad_) * Δxad_
        dtCFe[end] += neumann_right(DFeMg,CMg,ix,Δxad_) * Δxad_ +
                      neumann_right(DFeFe,CFe,ix,Δxad_) * Δxad_ +
                      neumann_right(DFeMn,CMn,ix,Δxad_) * Δxad_
        dtCMn[end] += neumann_right(DMnMg,CMg,ix,Δxad_) * Δxad_ +
                      neumann_right(DMnFe,CFe,ix,Δxad_) * Δxad_ +
                      neumann_right(DMnMn,CMn,ix,Δxad_) * Δxad_
    end
end

function semi_discretisation_diffusion_1D(du,u,p,t)

    @unpack D, D0, D_charact, Δxad_, bc_neumann = p.domain

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    dtCMg = @view du[:,1]
    dtCFe = @view du[:,2]
    dtCMn = @view du[:,3]

    # update diffusive parameters
    Diffusion_coef_1D!(D, CMg, CFe ,CMn, D0, D_charact)

    # semi-discretization
    stencil_diffusion_1D!(dtCMg, dtCFe, dtCMn, CMg, CFe ,CMn, D, Δxad_, bc_neumann)
end

