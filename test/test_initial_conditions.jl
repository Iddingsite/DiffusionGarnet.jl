
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
    C0 = ones(5) .* 0.1
    Lx = 10.0u"µm"
    tfinal = 1.0u"Myr"

    Grt_REE = Garnet.Grt_REE_Bloch2020_slow
    D = SetChemicalDiffusion(Grt_REE)

    IC1D = IC1DTrace(;C0, D, Lx, tfinal)

    @test IC1D.C0 == C0
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
    ICSph = ICSphTrace(;C0, D, Lr, tfinal)

    @test ICSph.C0 == C0
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

    IC2D = IC2DTrace(;C0=CMg0, D=D, Lx, Ly, tfinal, grt_boundary)

    @test IC2D.C0 == CMg0
    @test IC2D.Lx == ustrip(u"µm",Lx)
    @test IC2D.Ly == ustrip(u"µm",Ly)
    @test IC2D.tfinal == ustrip(u"Myr",tfinal)
    @test IC2D.D == D

    domain2D = Domain(IC2D, T, P)

    @test domain2D.L_charact == ustrip(u"µm",Lx)
    @test domain2D.t_charact == 1.0
    @test domain2D.tfinal_ad == 1.0
    @test domain2D.D_charact == 100.0

    C0 = ones(5, 5, 5) .* 0.1

    Lz = 10.0u"µm"
    tfinal = 1.0u"Myr"

    IC3D = IC3DTrace(;C0, D, Lx, Ly, Lz, tfinal)

    @test IC3D.C0 == C0
    @test IC3D.Lx == ustrip(u"µm",Lx)
    @test IC3D.Ly == ustrip(u"µm",Ly)
    @test IC3D.Lz == ustrip(u"µm",Lz)
    @test IC3D.tfinal == ustrip(u"Myr",tfinal)
    @test IC3D.D == D
    @test size(IC3D.C0, 1) == 5

    domain3D = Domain(IC3D, T, P)

    @test domain3D.L_charact == ustrip(u"µm",Lx)
    @test domain3D.t_charact == 1.0
    @test domain3D.tfinal_ad == 1.0
    @test domain3D.D_charact == 100.0
    @test domain3D.Δzad_ == 4.0
end
