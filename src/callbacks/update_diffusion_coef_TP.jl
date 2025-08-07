
"""
    update_diffusion_coef(integrator)

Callback function to update the diffusion coefficients at a given time from a new pressure and temperature. To use with the callback `PresetTimeCallback` (https://docs.sciml.ai/stable/basics/callbacks/#PresetTimeCallback-1).

Follows the syntax of callback functions defined by DiffEqCallbacks.jl (https://docs.sciml.ai/DiffEqCallbacks/stable/).
"""
function update_diffusion_coef(integrator)

    @unpack D0, P, T, fugacity_O2, time_update_ad, t_charact, diffcoef, D0_data = integrator.p.domain

    # find the index of the time_update_ad that is equal to t
    index = findfirst(x -> x == integrator.t, time_update_ad)

    nd = ndims(integrator.u)
    idx = ntuple(_ -> Colon(), nd - 1)
    CMg = @view integrator.u[idx..., 1]
    CFe = @view integrator.u[idx..., 2]
    CMn = @view integrator.u[idx..., 3]

    # update diffusion coefficients
    if index !== nothing

        for I in CartesianIndices(CMg)
            I_tuple = Tuple(I)  # convert CartesianIndex to Tuple for indexing
            # D0 is indexed as D0[:, I...]
            D0_view = @view D0[:, I_tuple...]

            # Extract local values at spatial index I
            cMg = CMg[I_tuple...]
            cFe = CFe[I_tuple...]
            cMn = CMn[I_tuple...]
            T_local = T[index]
            P_local = P[index]
            fO2_local = fugacity_O2[index]

            # Update diffusion coefficients at this point
            D_update!(D0_view, T_local, P_local, diffcoef, cMg, cFe, cMn, D0_data, fO2_local)
        end

        if integrator.t ≠ 0.0
            @info "New temperature and pressure: $(T[index]) °C and $(P[index]) kbar, updated at $(round((integrator.t * t_charact), digits=2)) Myr."
        end
    end
end
