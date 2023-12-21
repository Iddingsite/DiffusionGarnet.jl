using DiffusionGarnet
using Test
using DelimitedFiles

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
    @test IC2D.tfinal == ustrip(u"Myr",tfinal)

    CMg = ones(5, 5, 5) .* 0.1
    CFe = ones(5, 5, 5) .* 0.1
    CMn = ones(5, 5, 5) .* 0.1
    Lz = 10.0u"µm"

    IC3D = InitialConditions3D(CMg, CFe, CMn, Lx, Ly, Lz, tfinal)

    @test IC3D.CMg0 == CMg
    @test IC3D.Lz == ustrip(u"µm",Lz)
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

    using LinearAlgebra: norm


    data = DelimitedFiles.readdlm("./Data/1D/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

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

    sol = simulate(domain1D; progressbar=false)

    @test norm(sum.(sol[end][:,1] .+ sol[end][:,2] .+ sol[end][:,3])) ≈ 28.64886878627501
end

@testset "Spherical diffusion" begin

    using LinearAlgebra: norm


    data = DelimitedFiles.readdlm("./Data/1D/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

    Mg0 = reverse(data[1:size(data,1)÷2, 4])
    Fe0 = reverse(data[1:size(data,1)÷2, 2])
    Mn0 = reverse(data[1:size(data,1)÷2, 3])
    Ca0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lr = (data[end,1] - data[1,1])u"µm"
    tfinal = 15u"Myr"

    ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lr, tfinal)

    T = 900u"°C"
    P = 0.6u"GPa"

    domainSph = Domain(ICSph, T, P)

    sol = simulate(domainSph; progressbar=false)

    @test norm(sum.(sol[end][:,1] .+ sol[end][:,2] .+ sol[end][:,3])) ≈ 20.268803083443927
end

@testset "2D Diffusion" begin

    using LinearAlgebra: norm

    CMg = DelimitedFiles.readdlm("./Data/2D/Xprp_LR.txt", '\t', '\n', header=false)
    CFe = DelimitedFiles.readdlm("./Data/2D/Xalm_LR.txt", '\t', '\n', header=false)
    CMn = DelimitedFiles.readdlm("./Data/2D/Xsps_LR.txt", '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm("./Data/2D/Contour_LR.txt", '\t', '\n', header=false)

    Lx = 900.0u"µm"
    Ly = 900.0u"µm"
    tfinal = 1.0u"Myr"
    T = 900u"°C"
    P = 0.6u"GPa"

    IC2D = InitialConditions2D(CMg, CFe, CMn, Lx, Ly, tfinal; grt_boundary = grt_boundary)
    domain2D = Domain(IC2D, T, P)

    sol = simulate(domain2D; progressbar=false)

    @test norm(sol[end][:,:,1]) ≈ 12.783357041653609

end


@testset "Callback update D0" begin

    data = DelimitedFiles.readdlm("./Data/1D/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

    Mg0 = reverse(data[1:size(data,1)÷2, 4])
    Fe0 = reverse(data[1:size(data,1)÷2, 2])
    Mn0 = reverse(data[1:size(data,1)÷2, 3])
    Ca0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lx = Lr = (data[end,1] - data[1,1])u"µm"
    tfinal = 3u"Myr"

    ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lr, tfinal)
    IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)

    time_update = [0u"Myr", 2u"Myr"]
    T = [850u"°C", 600u"°C"]
    P = [0.5u"GPa", 0.3u"GPa"]

    domainSph = Domain(ICSph, T, P, time_update)
    domain1D = Domain(IC1D, T, P, time_update)

    @test domainSph.D0[1] ≈ 151880.41527919917
    @test domain1D.D0[1] ≈ 151880.41527919917

    CMg = DelimitedFiles.readdlm("./Data/2D/Xprp_LR.txt", '\t', '\n', header=false)
    CFe = DelimitedFiles.readdlm("./Data/2D/Xalm_LR.txt", '\t', '\n', header=false)
    CMn = DelimitedFiles.readdlm("./Data/2D/Xsps_LR.txt", '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm("./Data/2D/Contour_LR.txt", '\t', '\n', header=false)

    Lx = 900.0u"µm"
    Ly = 900.0u"µm"

    IC2D = InitialConditions2D(CMg, CFe, CMn, Lx, Ly, tfinal; grt_boundary = grt_boundary)

    domain2D = Domain(IC2D, T, P, time_update)

    @test domain2D.D0[1] ≈ 151880.41527919917


    @unpack time_update_ad = domainSph

    update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

    sol_sph = simulate(domainSph; callback=update_diffusion_coef_call, progressbar=false)
    sol_1D = simulate(domain1D; callback=update_diffusion_coef_call, progressbar=false)

    @unpack time_update_ad = domain2D
    update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

    sol_2D = simulate(domain2D; callback=update_diffusion_coef_call, progressbar=true)

    T=600  # in °C
    P=3  # in kbar
    D0 = zeros(Float64, 4)
    D_ini!(D0,T,P)

    @test sol_1D.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_sph.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_2D.prob.p.domain.D0[1] ≈ D0[1]
end

@testset "Callback output" begin

    using HDF5

    data_1D = DelimitedFiles.readdlm("./Data/1D/Data_Grt_1D.txt", '\t', '\n', header=true)[1]

    Mg0 = data_1D[:, 4]
    Fe0 = data_1D[:, 2]
    Mn0 = data_1D[:, 3]
    Ca0 = data_1D[:, 5]
    distance = data_1D[:, 1]
    Lx = (data_1D[end,1] - data_1D[1,1])u"µm"
    tfinal = 15u"Myr"

    IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)
    ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lx, tfinal)

    T = 700u"°C"
    P = 0.6u"GPa"

    domainSph = Domain(ICSph, T, P)
    domain1D = Domain(IC1D, T, P)

    Mg_2D = DelimitedFiles.readdlm("./Data/2D/Xprp_LR.txt", '\t', '\n', header=false)
    Fe_2D = DelimitedFiles.readdlm("./Data/2D/Xalm_LR.txt", '\t', '\n', header=false)
    Mn_2D = DelimitedFiles.readdlm("./Data/2D/Xsps_LR.txt", '\t', '\n', header=false)
    Lx = 900.0u"µm"
    Ly = 900.0u"µm"

    IC2D = InitialConditions2D(Mg_2D, Fe_2D, Mn_2D, Lx, Ly, tfinal)
    domain2D = Domain(IC2D, T, P)

    time_save = [0u"Myr", 2u"Myr", 5u"Myr", 15u"Myr"]

    save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain1D.t_charact, save_data)

    sol_1D = simulate(domain1D; callback=save_data_callback, path_save=(@__DIR__) * "/Grt_1D.h5", progressbar=false)

    sol_sph = simulate(domainSph; callback=save_data_callback, path_save=(@__DIR__) * "/Grt_Sph.h5", progressbar=false)

    save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain2D.t_charact, save_data)

    sol_2D = simulate(domain2D; callback=save_data_callback, path_save=(@__DIR__) * "/Grt_2D.h5", progressbar=false)

    h5open("./Grt_1D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == IC1D.CMg0
        @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == IC1D.CFe0
        @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == IC1D.CMn0
        @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == 1 .- IC1D.CMg0 .- IC1D.CFe0 .- IC1D.CMn0
        @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == sol_1D[end][:,1]
        @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == sol_1D[end][:,2]
        @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == sol_1D[end][:,3]
    end

    h5open("./Grt_Sph.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == ICSph.CMg0
        @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == ICSph.CFe0
        @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == ICSph.CMn0
        @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == 1 .- ICSph.CMg0 .- ICSph.CFe0 .- ICSph.CMn0
        @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == sol_sph[end][:,1]
        @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == sol_sph[end][:,2]
        @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == sol_sph[end][:,3]
    end

    h5open("./Grt_2D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"])' == IC2D.CMg0
        @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"])' == IC2D.CFe0
        @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"])' == IC2D.CMn0
        @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"])' == 1 .- IC2D.CMg0 .- IC2D.CFe0 .- IC2D.CMn0
        @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"])' == sol_2D[end][:,:,1]
        @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"])' == sol_2D[end][:,:,2]
        @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"])' == sol_2D[end][:,:,3]
    end

    # delete files "Grt_1D.h5" and "Grt_Sph.h5" if it exists
    if isfile("./Grt_1D.h5")
        rm("./Grt_1D.h5")
    end

    if isfile("./Grt_Sph.h5")
        rm("./Grt_Sph.h5")
    end

    if isfile("./Grt_2D.h5")
        rm("./Grt_2D.h5")
    end

end