function D_update!(D0,T,P)
    R = 8.314462618  # constant of gas in J mol−1 K−1

    # Magnesium
    D0Mg = 1.1 * 1e-3 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Mg = 67997 * 4.1855 # activation energy at 1 bar in J / mol
    ΔV⁺Mg = 5.3  # activation volume in cm3 / mol

    # Iron
    D0Fe = 6.4 * 1e-4 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Fe = 65824 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Fe = 5.6  # activation volume in cm3 / mol

    # Manganese
    D0Mn = 5.1 * 1e-4 * 1e8  # pre-exponential constant in µm2 / s
    Eₐ_Mn = 60569 * 4.1855  # activation energy at 1 bar in J / mol
    ΔV⁺Mn = 6.0  # activation volume in cm3 / mol

    DMg = D0Mg * exp(- (Eₐ_Mg + (100 * (P-0.001) * ΔV⁺Mg)) / (R * (T+273.15)))  # in µm2 / s
    DFe = D0Fe * exp(- (Eₐ_Fe + (100 * (P-0.001) * ΔV⁺Fe)) / (R * (T+273.15)))  # in µm2 / s
    DMn = D0Mn * exp(- (Eₐ_Mn + (100 * (P-0.001) * ΔV⁺Mn)) / (R * (T+273.15)))  # in µm2 / s
    DCa = 0.5 * DFe


    D0 .= [DMg, DFe, DMn, DCa] .* (365.25 * 24 * 3600 * 1e6)  # in years
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
            println("New temperature and pressure: $(T[index]) °C and $(P[index]) kbar, updated at $(round((integrator.t * t_charact), digits=2)) Myr.")
        end
    end
end
