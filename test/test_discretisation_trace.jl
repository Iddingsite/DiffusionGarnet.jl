
@testset "Trace element diffusion 1D and Spherical" begin

    nx = 400
    L = 0.8u"mm"
    Grt_Mg = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_Mg)

    # Gaussian initial condition
    C = zeros(nx)
    C[1:nx÷2] .= 4.5
    C[nx÷2+1:end] .= 1.2

    # 1D case
    IC1D = IC1DTrace(;C0=C, D=D, Lx=L, tfinal=1u"Myr")
    domain1D = Domain(IC1D, 800u"°C", 1u"GPa")
    sol1D = simulate(domain1D; progress=false, abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=true)
    @test sol1D.u[1] == C
    @test length(sol1D.u[end]) == nx
    # use norm 
    @test norm(sol1D.u[end]) ≈ 60.67095414913916

    # Spherical case
    r = range(0u"mm", length=nx, stop=L)
    ICsph = ICSphTrace(;C0=C, D=D, Lr=L, r=r, tfinal=1u"Myr")
    domainsph = Domain(ICsph, 800u"°C", 1u"GPa")
    solsph = simulate(domainsph; progress=false, abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=true)
    @test solsph.u[1] == C
    @test length(solsph.u[end]) == nx
    @test norm(solsph.u[end]) ≈ 43.00182957068845
end

@testset "Trace element diffusion 2D" begin

    root = @__DIR__
    path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
    path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

    CMg0 = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

    # multiply by 1000 to simulate ppm
    CMg0 .= CMg0 .* 1e5

    Lx = 10.0u"µm"
    Ly = 10.0u"µm"
    tfinal = 1.0u"Myr"
    Grt_Mg = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_Mg)

    IC2D = IC2DTrace(;C0=CMg0, D=D, Lx, Ly, tfinal, grt_boundary)
    domain2D = Domain(IC2D, 800u"°C", 1u"GPa")

    sol2D = simulate(domain2D; progress=false, abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=true)

    @test sol2D.u[1] == CMg0
    @test size(sol2D.u[end]) == size(CMg0)#
    @test norm(sol2D.u[end]) ≈ 1.2784768016330013e6
end

@testset "Trace element diffusion 3D" begin


    # use JLD2 to load data
    root = @__DIR__
    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")

    file = jldopen(path_3D, "r")
    Mg0 = file["Mg0"]
    grt_boundary = file["grt_boundary"]
    close(file)

    Grt_Mg = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_Mg)

    Lx = 9000.0u"µm"
    Ly = 9000.0u"µm"
    Lz = 9000.0u"µm"
    tfinal = 0.2u"Myr"
    T = 900u"°C"
    P = 0.6u"GPa"

    Mg0 .= Mg0 .* 1e5

    IC3D = IC3DTrace(;C0=Mg0, D=D, Lx, Ly, Lz, tfinal, grt_boundary)
    domain3D = Domain(IC3D, T, P)

    sol3D = simulate(domain3D; progress=false, abstol=1e-6, reltol=1e-6, save_everystep=false, save_start=true);

    @test sol3D.u[1] == Mg0
    @test size(sol3D.u[end]) == size(Mg0)
    @test norm(sol3D.u[end]) ≈ 3.675704525002046e7
end
