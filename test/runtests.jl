using DiffusionGarnet
using Test

@testset "initial conditions" begin
    # 1D, test geometry.jl
    CMg = ones(5) .* 0.1
    CFe = ones(5) .* 0.1
    CMn = ones(5) .* 0.1
    Lx = 10.0u"µm"
    tfinal = 1.0u"Myr"

    IC1D = InitialConditions1D(CMg, CFe, CMn, Lx, tfinal)

    @test IC1D.CMg0 == CMg
    @test IC1D.Lx == ustrip(u"µm",Lx)
    @test IC1D.tfinal == ustrip(u"Myr",tfinal)

    ICSph = InitialConditionsSpherical(CMg, CFe, CMn, Lx, tfinal)

    @test ICSph.CMg0 == CMg
    @test ICSph.Lr == ustrip(u"µm",Lx)
    @test ICSph.tfinal == ustrip(u"Myr",tfinal)

    CMg = ones(5, 5) .* 0.1
    CFe = ones(5, 5) .* 0.1
    CMn = ones(5, 5) .* 0.1
    Ly = 10.0u"µm"

    IC2D = InitialConditions2D(CMg, CFe, CMn, Lx, Ly, tfinal)

    @test IC2D.CMg0 == CMg
    @test IC2D.Ly == ustrip(u"µm",Ly)
    @test IC2D.grid.x == IC2D.x' .* ones(size(CMg,2))
    @test IC2D.tfinal == ustrip(u"Myr",tfinal)

    CMg = ones(5, 5, 5) .* 0.1
    CFe = ones(5, 5, 5) .* 0.1
    CMn = ones(5, 5, 5) .* 0.1
    Lz = 10.0u"µm"

    IC3D = InitialConditions3D(CMg, CFe, CMn, Lx, Ly, Lz, tfinal)

    @test IC3D.CMg0 == CMg
    @test IC3D.Lz == ustrip(u"µm",Lz)
    @test IC3D.grid.x == IC3D.x' .* ones(size(CMg,2), size(CMg,3))
    @test IC3D.tfinal == ustrip(u"Myr",tfinal)

    T = 650u"°C"
    P = 2u"GPa"
    domain1D = Domain(IC1D, T, P)
    @test domain1D.L_charact == ustrip(u"µm",Lx)
    @test domain1D.t_charact ≈ 0.22518558662307234
    @test domain1D.tfinal_ad ≈ 4.44078155709785
    @test domain1D.Δxad_ == 4.0

    domainSph = Domain(ICSph, T, P)
    @test domainSph.L_charact == ustrip(u"µm",Lx)
    @test domainSph.t_charact ≈ 0.22518558662307234
    @test domainSph.tfinal_ad ≈ 4.44078155709785
    @test domainSph.Δrad_ == 4.0
    @test domainSph.r_ad[end] == 1.0

    domain2D = Domain(IC2D, T, P)
    @test domain2D.L_charact == ustrip(u"µm",Lx)
    @test domain2D.Δyad_ == 4.0

    domain3D = Domain(IC3D, T, P)
    @test domain3D.L_charact == ustrip(u"µm",Lx)
    @test domain3D.Δzad_ == 4.0
end

@testset "diffusion coefficients" begin
    # test Domain
    T=650  # in °C
    P=2  # in kbar
    D0 = zeros(Float64, 4)
    D_ini!(D0,T,P)

    @test D0[1] ≈ 2.383676419323230e2
    @test D0[2] ≈ 4.500484945403809e+02
    @test D0[3] ≈ 6.232092221232668e+03
    @test D0[4] ≈ 2.250242472701905e+02
end

@testset "1D diffusion" begin

    data = DelimitedFiles.readdlm("./Data_Grt_1D.txt", '\t', '\n', header=true)[1]

    Mg0 = data[:, 4]
    Fe0 = data[:, 2]
    Mn0 = data[:, 3]
    Ca0 = data[:, 5]
    distance = data[:, 1]
    Lx = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domain1D = Domain(IC1D, T, P)

    sol = simulate(domain1D)

    @test norm(sum.(sol[end][:,1] .+ sol[end][:,2] .+ sol[end][:,3])) ≈ 28.64886878627501
end

@testset "Spherical diffusion" begin

    data = DelimitedFiles.readdlm("./Data_Grt_1D.txt", '\t', '\n', header=true)[1]

    Mg0 = reverse(data[1:size(data,1)÷2, 4])
    Fe0 = reverse(data[1:size(data,1)÷2, 2])
    Mn0 = reverse(data[1:size(data,1)÷2, 3])
    Ca0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lx = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lx, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domainSph = Domain(ICSph, T, P)

    sol = simulate(domainSph)

    @test norm(sum.(sol[end][:,1] .+ sol[end][:,2] .+ sol[end][:,3])) ≈ 20.268803083443927
end
