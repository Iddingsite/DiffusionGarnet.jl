function D_update!(D0,T,P, fugacity_O2=1e-25)

    Grt_Mg = Grt_Mg_Chakraborty1992
    Grt_Fe = Grt_Fe_Chakraborty1992
    Grt_Mn = Grt_Mn_Chakraborty1992

    Grt_Mg = SetChemicalDiffusion(Grt_Mg)
    Grt_Fe = SetChemicalDiffusion(Grt_Fe)
    Grt_Mn = SetChemicalDiffusion(Grt_Mn)

    T_K = (T+273.15) .* 1u"K"
    P_kbar = P * 1u"kbar"

    DMg = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Mg, T = T_K, P = P_kbar)))
    DFe = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Fe, T = T_K, P = P_kbar)))
    DMn = ustrip(uconvert(u"µm^2/Myr",compute_D(Grt_Mn, T = T_K, P = P_kbar)))
    DCa = 0.5 * DFe

    fugacity_ratio = fugacity_O2/1e-25  # current fO2 over fO2 buffered with graphite

    # after Chakraborty and Ganguly, 1991 (page 142, equation 2)
    DMg = exp(log(DMg) + 1/6 * log(fugacity_ratio))
    DFe = exp(log(DFe) + 1/6 * log(fugacity_ratio))
    DMn = exp(log(DMn) + 1/6 * log(fugacity_ratio))
    DCa = exp(log(DCa) + 1/6 * log(fugacity_ratio))

    D0 .= (DMg, DFe, DMn, DCa)   # in µm^2/Myr
end

"""
    update_diffusion_coef(integrator)

Callback function to update the diffusion coefficients at a given time from a new pressure and temperature. To use with the callback `PresetTimeCallback` (https://docs.sciml.ai/stable/basics/callbacks/#PresetTimeCallback-1).

Follows the syntax of callback functions defined by DiffEqCallbacks.jl (https://docs.sciml.ai/DiffEqCallbacks/stable/).
"""
function update_diffusion_coef(integrator)

    @unpack D0, P, T, time_update_ad, t_charact = integrator.p.domain

    # find the index of the time_update_ad that is equal to t
    index = findfirst(x -> x == integrator.t, time_update_ad)

    #! integrator.opts.callback

    # update diffusion coefficients
    if index !== nothing
        D_update!(D0, T[index], P[index])

        if integrator.t ≠ 0.0
            @info "New temperature and pressure: $(T[index]) °C and $(P[index]) kbar, updated at $(round((integrator.t * t_charact), digits=2)) Myr."
        end
    end
end
