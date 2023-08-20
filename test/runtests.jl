using DiffusionGarnet
using Test

@testset "initial conditions" begin
    # 1D, test geometry.jl
    CMg = ones(5)
    CFe = ones(5)
    CMn = ones(5)
    Lx = 10.0u"km"
    tfinal = 10.0u"s"

    Ini1D = InitialConditions1D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, tfinal=tfinal)

    @test Ini1D.CMg0 == CMg
    @test Ini1D.Lx == ustrip(u"m",Lx)
    @test Ini1D.tfinal == ustrip(u"s",tfinal)

    CMg = ones(5, 5)
    CFe = ones(5, 5)
    CMn = ones(5, 5)
    Ly = 10.0u"km"

    Ini2D = InitialConditions2D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, Ly=Ly, tfinal=tfinal)

    @test Ini2D.CMg0 == CMg
    @test Ini2D.Ly == ustrip(u"m",Ly)
    @test Ini2D.grid.x == Ini2D.x' .* ones(size(CMg,2))
    @test Ini2D.tfinal == ustrip(u"s",tfinal)

    CMg = ones(5, 5, 5)
    CFe = ones(5, 5, 5)
    CMn = ones(5, 5, 5)
    Lz = 10.0u"km"

    Ini3D = InitialConditions3D(CMg0=CMg, CFe0=CFe, CMn0=CMn, Lx=Lx, Ly=Ly, Lz=Lz, tfinal=tfinal)

    @test Ini3D.CMg0 == CMg
    @test Ini3D.Lz == ustrip(u"m",Lz)
    @test Ini3D.grid.x == Ini3D.x' .* ones(size(CMg,2), size(CMg,3))
end
