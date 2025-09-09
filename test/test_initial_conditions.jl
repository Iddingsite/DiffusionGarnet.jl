
@testset "All tests" begin

    @testset "initial conditions Major" begin

        CMg = ones(5) .* 0.1
        CFe = ones(5) .* 0.1
        CMn = ones(5) .* 0.1
        Lx = 10.0u"µm"
        tfinal = 1.0u"Myr"

        IC1D = IC1DMajor(;CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx, tfinal)

        @test IC1D.CMg0 == CMg
        @test IC1D.Lx == ustrip(u"µm",Lx)
        @test IC1D.tfinal == ustrip(u"Myr",tfinal)

        ICSph = ICSphMajor(;CMg0=CMg, CFe0=CFe, CMn0=CMn, Lr=Lx, tfinal)

        @test ICSph.CMg0 == CMg
        @test ICSph.Lr == ustrip(u"µm",Lx)
        @test ICSph.tfinal == ustrip(u"Myr",tfinal)

        CMg0 = ones(5, 5) .* 0.1
        CFe0 = ones(5, 5) .* 0.1
        CMn0 = ones(5, 5) .* 0.1
        Ly = 10.0u"µm"

        IC2D = IC2DMajor(;CMg0, CFe0, CMn0, Lx, Ly, tfinal)

        @test IC2D.CMg0 == CMg0
        @test IC2D.Ly == ustrip(u"µm",Ly)
        @test IC2D.tfinal == ustrip(u"Myr",tfinal)

        CMg0 = ones(5, 5, 5) .* 0.1
        CFe0 = ones(5, 5, 5) .* 0.1
        CMn0 = ones(5, 5, 5) .* 0.1
        Lz = 10.0u"µm"

        IC3D = IC3DMajor(;CMg0, CFe0, CMn0, Lx, Ly, Lz, tfinal)

        @test IC3D.CMg0 == CMg0
        @test IC3D.Lz == ustrip(u"µm",Lz)
        @test IC3D.tfinal == ustrip(u"Myr",tfinal)

        T = 650u"°C"
        P = 2u"GPa"
        domain1D = Domain(IC1D, T, P)
        @test domain1D.L_charact == ustrip(u"µm",Lx)
        @test domain1D.t_charact ≈ 1
        @test domain1D.tfinal_ad ≈ 1
        @test domain1D.Δxad_[1] == 4.0

        domainSph = Domain(ICSph, T, P)
        @test domainSph.L_charact == ustrip(u"µm",Lx)
        @test domainSph.t_charact ≈ 1
        @test domainSph.tfinal_ad ≈ 1
        @test domainSph.Δr_ad_[1] == 4.0
        @test domainSph.r_ad[end] == 1.0

        domain2D = Domain(IC2D, T, P)
        @test domain2D.L_charact == ustrip(u"µm",Lx)
        @test domain2D.Δyad_ == 4.0

        domain3D = Domain(IC3D, T, P)
        @test domain3D.L_charact == ustrip(u"µm",Lx)
        @test domain3D.Δzad_ == 4.0
    end


    @testset "initial conditions Trace" begin
        C = ones(5) .* 0.1
        Lx = 10.0u"µm"
        tfinal = 1.0u"Myr"

        Grt_REE = Garnet.Grt_REE_Bloch2020_slow
        D = SetChemicalDiffusion(Grt_REE)

        IC1D = IC1DTrace(;C, D, Lx, tfinal)

        @test IC1D.C == C
        @test IC1D.Lx == ustrip(u"µm",Lx)
        @test IC1D.tfinal == ustrip(u"Myr",tfinal)
        @test IC1D.D == D

        T = 650u"°C"
        P = 2u"GPa"
        domain1D = Domain(IC1D, T, P)
        @test domain1D.L_charact == ustrip(u"µm",Lx)
        @test domain1D.t_charact == 1.0
        @test domain1D.tfinal_ad == 1.0
        @test domain1D.D_charact == 100.0
        @test domain1D.Δxad_[1] == 4.0

        Lr = Lx
        ICSph = ICSphTrace(;C, D, Lr, tfinal)

        @test ICSph.C == C
        @test ICSph.Lr == ustrip(u"µm",Lr)
        @test ICSph.tfinal == ustrip(u"Myr",tfinal)
        @test ICSph.D == D

        domainSph = Domain(ICSph, T, P)
        @test domainSph.L_charact == ustrip(u"µm",Lr)
        @test domainSph.t_charact == 1.0
        @test domainSph.tfinal_ad == 1.0
        @test domainSph.D_charact == 100.0
        @test domainSph.Δr_ad_[1] == 4.0

        root = @__DIR__
        path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
        path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

        CMg0 = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
        grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

        # multiply by 1000 to simulate ppm
        CMg0 .= CMg0 .* 1e5

        Ly = 10.0u"µm"
        tfinal = 1.0u"Myr"

        IC2D = IC2DTrace(;C=CMg0, D=D, Lx, Ly, tfinal, grt_boundary)

        @test IC2D.C == CMg0
        @test IC2D.Lx == ustrip(u"µm",Lx)
        @test IC2D.Ly == ustrip(u"µm",Ly)
        @test IC2D.tfinal == ustrip(u"Myr",tfinal)
        @test IC2D.D == D

        domain2D = Domain(IC2D, T, P)

        @test domain2D.L_charact == ustrip(u"µm",Lx)
        @test domain2D.t_charact == 1.0
        @test domain2D.tfinal_ad == 1.0
        @test domain2D.D_charact == 100.0

        C = ones(5, 5, 5) .* 0.1

        Lz = 10.0u"µm"
        tfinal = 1.0u"Myr"

        IC3D = IC3DTrace(;C, D, Lx, Ly, Lz, tfinal)

        @test IC3D.C == C
        @test IC3D.Lx == ustrip(u"µm",Lx)
        @test IC3D.Ly == ustrip(u"µm",Ly)
        @test IC3D.Lz == ustrip(u"µm",Lz)
        @test IC3D.tfinal == ustrip(u"Myr",tfinal)
        @test IC3D.D == D
        @test size(IC3D.C, 1) == 5

        domain3D = Domain(IC3D, T, P)

        @test domain3D.L_charact == ustrip(u"µm",Lx)
        @test domain3D.t_charact == 1.0
        @test domain3D.tfinal_ad == 1.0
        @test domain3D.D_charact == 100.0
        @test domain3D.Δzad_ == 4.0
    end

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


        DiffusionGarnet.D_update!(D0,T_K,P_kbar,D0_data)

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
end






