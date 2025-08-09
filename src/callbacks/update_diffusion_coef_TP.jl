
"""
    D_update!(D0, T_K, P_kbar, diffcoef, CMg, CFe, CMn, X, fugacity_O2=1e-25NoUnits)

Update the diffusion coefficients `D0` based on the temperature `T_K`, pressure `P_kbar`, and the chemical compositions `CMg`, `CFe`, and `CMn`. The optional parameter `fugacity_O2` is used to set the oxygen fugacity, defaulting to 1e-25 Pa (graphite buffer).
"""
function D_update!(D0, T_K, P_kbar, D0_data, fugacity_O2=1e-25NoUnits)  # by defaut 1e-25 Pa is graphite buffer

    DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2)))
    DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2)))
    DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(D0_data.Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2)))

    DCa = 0.5 * DFe

    D0 .= (DMg, DFe, DMn, DCa)   # in µm^2/Myr
end


function update_diffusion_coef(integrator)

    @unpack D0, P, T, fugacity_O2, time_update_ad, D0_data, diffcoef = integrator.p.domain

    # find the index of the time_update_ad that is equal to t
    index = findfirst(x -> x == integrator.t, time_update_ad)

    #! integrator.opts.callback

    # update diffusion coefficients
    if index !== nothing

        # only for Chakraborty1992
        if diffcoef == 1

            T_K = (T[index]+273.15) * 1u"K"
            P_kbar = P[index] * 1u"kbar"
            fO2 = (fugacity_O2[index])NoUnits  # default value for graphite

            D_update!(D0, T_K, P_kbar, D0_data, fO2)
        end

        if integrator.t ≠ 0.0
            @info "New temperature and pressure: $(round(T[index],digits=2)) °C and $(round(P[index]*0.1, digits=2)) GPa, updated at $(round((integrator.t), digits=2)) Myr."
        end
    end
end
