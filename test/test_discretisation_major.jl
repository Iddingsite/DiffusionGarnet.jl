
@testset "1D diffusion" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    CMg0 = data[:, 4]
    CFe0 = data[:, 2]
    CMn0 = data[:, 3]
    CCa0 = data[:, 5]
    distance = data[:, 1]
    Lx = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    IC1D = IC1DMajor(;CMg0, CFe0, CMn0, Lx, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domain1D = Domain(IC1D, T, P)

    sol = simulate(domain1D; progress=false, save_everystep=false, save_start=false)

    @test norm(sum.(sol.u[end][:,1] .+ sol.u[end][:,2] .+ sol.u[end][:,3])) ≈ 28.653159661729227 rtol=1e-5
end

@testset "Spherical diffusion" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    Mg0 = reverse(data[1:size(data,1)÷2, 4])
    Fe0 = reverse(data[1:size(data,1)÷2, 2])
    Mn0 = reverse(data[1:size(data,1)÷2, 3])
    Ca0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lr = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    ICSph = ICSphMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lr, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domainSph = Domain(ICSph, T, P)

    sol = simulate(domainSph; progress=false, save_everystep=false, save_start=false)

    @test norm(sum.(sol.u[end][:,1] .+ sol.u[end][:,2] .+ sol.u[end][:,3])) ≈ 20.272049461615836 rtol=1e-5
end

@testset "2D Diffusion" begin

    root = @__DIR__
    path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
    path_2D_Fe = joinpath(root, "Data", "2D", "Xalm_LR.txt")
    path_2D_Mn = joinpath(root, "Data", "2D", "Xsps_LR.txt")
    path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

    CMg0 = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
    CFe0 = DelimitedFiles.readdlm(path_2D_Fe, '\t', '\n', header=false)
    CMn0 = DelimitedFiles.readdlm(path_2D_Mn, '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

    Lx = 900.0u"µm"
    Ly = 900.0u"µm"
    tfinal = 2.0u"Myr"
    T = 900u"°C"
    P = 0.6u"GPa"

    IC2D = IC2DMajor(;CMg0, CFe0, CMn0, Lx, Ly, tfinal, grt_boundary)
    domain2D = Domain(IC2D, T, P)

    sol = simulate(domain2D; save_everystep=false, save_start=false)

    @test norm(sol.u[end][:,:,1]) ≈ 12.783357393401664 rtol=1e-5
end

@testset "3D Diffusion" begin

    # use JLD2 to load data
    root = @__DIR__
    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")

    file = jldopen(path_3D, "r")
    @unpack Mg0, Fe0, Mn0, Ca0, grt_boundary = file
    close(file)

    Lx = 9000.0u"µm"
    Ly = 9000.0u"µm"
    Lz = 9000.0u"µm"
    tfinal = 0.2u"Myr"
    T = 900u"°C"
    P = 0.6u"GPa"

    IC3D = IC3DMajor(;CMg0=Mg0, CFe0=Fe0, CMn0=Mn0, Lx, Ly, Lz, tfinal, grt_boundary)
    domain3D = Domain(IC3D, T, P)

    sol = simulate(domain3D; save_everystep=false, save_start=false);

    @test norm(sol.u[end][:,:,:,1]) ≈ 371.1775756471261
end
