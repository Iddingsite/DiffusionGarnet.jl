
@testset "diffusion coefficients" begin
    # test Domain
    T=650  # in °C
    P=2  # in kbar
    diffcoef = 1  # Chakraborty and Ganguly 1992 diffusion coefficients
    CMg0 = 0.1
    CFe0 = 0.1
    CMn0 = 0.1

    Grt_Mg = SetChemicalDiffusion(Garnet.Grt_Mg_Chakraborty1992)
    Grt_Fe = SetChemicalDiffusion(Garnet.Grt_Fe_Chakraborty1992)
    Grt_Mn = SetChemicalDiffusion(Garnet.Grt_Mn_Chakraborty1992)

    D0_data = (Grt_Mg=Grt_Mg, Grt_Fe=Grt_Fe, Grt_Mn=Grt_Mn)

    D0 = zeros(Float64, 4)

    T_K = (T + 273.15) * u"K"
    P_kbar = P * u"kbar"
    fugacity_O2 = 1e-25NoUnits  # default value

    D0[1] = ustrip(u"µm^2/Myr", compute_D(Grt_Mg, T = T_K, P = P_kbar, fO2 = fugacity_O2))
    D0[2] = ustrip(u"µm^2/Myr", compute_D(Grt_Fe, T = T_K, P = P_kbar, fO2 = fugacity_O2))
    D0[3] = ustrip(u"µm^2/Myr", compute_D(Grt_Mn, T = T_K, P = P_kbar, fO2 = fugacity_O2))
    D0[4] = 0.5 * D0[2]  # DCa = 0.5 * DFe

    @test D0[1] ≈ 241.5563196143817
    @test D0[2] ≈ 455.8752224396371
    @test D0[3] ≈ 6306.29872433989
    @test D0[4] ≈ 227.93761121981854

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    Mg0 = data[:, 4]
    Fe0 = data[:, 2]
    Mn0 = data[:, 3]
    Ca0 = data[:, 5]
    distance = data[:, 1]
    Lx = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    IC1D = IC1DMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lx, tfinal)
    ICSph = ICSphMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lr=Lx, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domain1D = Domain(IC1D, T, P)
    domainSph = Domain(ICSph, T, P)

    @unpack D, D0, D_charact, Δxad_, bc_neumann = domain1D

    DiffusionGarnet.Diffusion_coef_1D_major!(domain1D.D, Mg0, Fe0, Mn0, domain1D.D0, domain1D.D_charact, domain1D, 0)
    DiffusionGarnet.Diffusion_coef_1D_major!(domainSph.D, Mg0, Fe0, Mn0, domainSph.D0, domainSph.D_charact, domainSph, 0)

    @test domain1D.D.DMgMg[1] ≈ 0.00688293300181906
    @test domain1D.D.DMgFe[1] ≈ -0.0006387171241768935
    @test domain1D.D.DMgMn[1] ≈ -0.008823562361791238
    @test domain1D.D.DFeFe[1] ≈ 0.005961179130686035
    @test domain1D.D.DFeMg[1] ≈ -0.0016564504861627976
    @test domain1D.D.DFeMn[1] ≈ -0.06043590833889515
    @test domain1D.D.DMnMn[1] ≈ 0.0747689622855202
    @test domain1D.D.DMnMg[1] ≈ -4.912582546362481e-5
    @test domain1D.D.DMnFe[1] ≈ -0.00012974512729293407

    @test domainSph.D.DMgMg[1] ≈ domain1D.D.DMgMg[1]
    @test domainSph.D.DMgFe[1] ≈ domain1D.D.DMgFe[1]
    @test domainSph.D.DMgMn[1] ≈ domain1D.D.DMgMn[1]
    @test domainSph.D.DFeFe[1] ≈ domain1D.D.DFeFe[1]
    @test domainSph.D.DFeMg[1] ≈ domain1D.D.DFeMg[1]
    @test domainSph.D.DFeMn[1] ≈ domain1D.D.DFeMn[1]
    @test domainSph.D.DMnMn[1] ≈ domain1D.D.DMnMn[1]
    @test domainSph.D.DMnMg[1] ≈ domain1D.D.DMnMg[1]
    @test domainSph.D.DMnFe[1] ≈ domain1D.D.DMnFe[1]

    # check diffusion coefficients from Bloch et al. 2020 from GeoParams
    Grt_Mg = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_Mg)

    D =  ustrip(u"m^2/s", GeoParams.compute_D(D, T = 800u"°C", P = 1u"GPa"))  # diffusion coefficient in µm^2/Myr
    D_paper = 10^(-10.24 - (221057) / (log(10) * 8.31446261815324 * (800+ 273.15)))

    @test D ≈ D_paper
end
