
"""
    D_update!(Domain, IC::T, T_K, P_kbar, fO2) where T <: InitialConditionsMajor

Update the diffusion coefficients `D0` based on the temperature `T_K`, pressure `P_kbar` and oxygen fugacity `fO2`
"""
function D_update!(Domain, IC::T, T_K, P_kbar, fO2) where T <: InitialConditionsMajor

    (; D0, D0_data, diffcoef) = Domain

    if diffcoef == 1 # because others are updated everytime step anyway
        DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fO2)))
        DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fO2)))
        DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fO2)))

        DCa = 0.5 * DFe

        D0 .= (DMg, DFe, DMn, DCa)   # in µm^2/Myr
    end

end

"""
    D_update!(Domain, IC::T, T_K, P_kbar, fO2) where T <: InitialConditionsTrace

Update the diffusion coefficients `D0` based on the temperature `T_K`, pressure `P_kbar` and oxygen fugacity `fO2`
"""
function D_update!(Domain, IC::T, T_K, P_kbar, fO2) where T <: InitialConditionsTrace

    (; fugacity_O2) = Domain

    Domain.D[1] = ustrip(u"µm^2/Myr", compute_D(IC.D, T = T_K, P = P_kbar, fugacity_O2 = fO2))
end

"""
    update_diffusion_coef(integrator)

Callback function to update the diffusion coefficients during the simulation based on the time steps defined in `time_update_ad`. Works for both major and trace elements.
"""
function update_diffusion_coef(integrator)

    (; P, T, fugacity_O2, time_update_ad, IC) = integrator.p.domain

    # find the index of the time_update_ad that is equal to t
    index = findfirst(x -> x == integrator.t, time_update_ad)

    # update diffusion coefficients
    if index !== nothing

        T_K = (T[index]+273.15) * 1u"K"
        P_kbar = P[index] * 1u"kbar"
        fO2 = (fugacity_O2[index])NoUnits  # default value for graphite

        D_update!(integrator.p.domain, IC, T_K, P_kbar, fO2)

        if integrator.t ≠ 0.0
            @info "New temperature and pressure: $(round(T[index],digits=2)) °C and $(round(P[index]*0.1, digits=2)) GPa, updated at $(round((integrator.t), digits=2)) Myr."
        end
    end
end
