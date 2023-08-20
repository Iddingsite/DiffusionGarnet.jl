using DiffusionGarnet
using Test

@testset "initial conditions" begin
    # 1D, test geometry.jl
    CMg = ones(5)
    CFe = ones(5)
    CMn = ones(5)
    Lx = 10.0u"km"
    Lx = 10.0u"m"
    tfinal = 1.0u"Myr"

    IC1D = InitialConditions1D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, tfinal=tfinal)

    @test IC1D.CMg0 == CMg
    @test IC1D.Lx == ustrip(u"m",Lx)
    @test IC1D.tfinal == ustrip(u"Myr",tfinal)

    CMg = ones(5, 5)
    CFe = ones(5, 5)
    CMn = ones(5, 5)
    Ly = 10.0u"km"

    IC2D = InitialConditions2D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, Ly=Ly, tfinal=tfinal)

    @test IC2D.CMg0 == CMg
    @test IC2D.Ly == ustrip(u"m",Ly)
    @test IC2D.grid.x == IC2D.x' .* ones(size(CMg,2))
    @test IC2D.tfinal == ustrip(u"Myr",tfinal)

    CMg = ones(5, 5, 5)
    CFe = ones(5, 5, 5)
    CMn = ones(5, 5, 5)
    Lz = 10.0u"km"

    IC3D = InitialConditions3D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, Ly=Ly, Lz=Lz, tfinal=tfinal)

    @test IC3D.CMg0 == CMg
    @test IC3D.Lz == ustrip(u"m",Lz)
    @test IC3D.grid.x == IC3D.x' .* ones(size(CMg,2), size(CMg,3))
    @test IC3D.tfinal == ustrip(u"Myr",tfinal)

    T = 650u"°C"
    P = 2u"GPa"
    domain1D = Domain1D(IC=IC1D, T=T, P=P)
    @test domain1D.L_charact == ustrip(u"m",Lx)
    @test domain1D.t_charact ≈ 0.22518558662307234
    @test domain1D.tfinal_ad ≈ 4.44078155709785
    @test domain1D.Δxad_ == 5.0

end

@testset "diffusion coefficients" begin
    # test Domain
    T=650
    P=2
    D0 = zeros(Float64, 4)
    D_ini!(D0,T,P)

    @test D0[1] ≈ 2.383676419323230e2
    @test D0[2] ≈ 4.500484945403809e+02
    @test D0[3] ≈ 6.232092221232668e+03
    @test D0[4] ≈ 2.250242472701905e+02
end
