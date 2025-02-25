using DiffusionGarnet
using Test
using DelimitedFiles
using JLD2
using HDF5
using LinearAlgebra: norm

@testset "All tests" begin

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
        @test domain1D.t_charact ≈ 0.2224999357930375
        @test domain1D.tfinal_ad ≈ 4.494383319418882
        @test domain1D.Δxad_ == 4.0

        domainSph = Domain(ICSph, T, P)
        @test domainSph.L_charact == ustrip(u"µm",Lx)
        @test domainSph.t_charact ≈ 0.2224999357930375
        @test domainSph.tfinal_ad ≈ 4.494383319418882
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
        DiffusionGarnet.D_ini!(D0,T,P)

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

        IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)
        ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lx, tfinal)

        T = 900u"°C"
        P = 0.6u"GPa"

        domain1D = Domain(IC1D, T, P)
        domainSph = Domain(ICSph, T, P)

        @unpack D, D0, D_charact, Δxad_, bc_neumann = domain1D

        DiffusionGarnet.Diffusion_coef_1D_major!(domain1D.D, Mg0, Fe0, Mn0, domain1D.D0, domain1D.D_charact)
        DiffusionGarnet.Diffusion_coef_spherical_major!(domainSph.D, Mg0, Fe0, Mn0, domainSph.D0, domainSph.D_charact)

        @test domain1D.D.DMgMg[1] ≈ 0.27756534824753715
        @test domain1D.D.DMgFe[1] ≈ -0.02575729575705166
        @test domain1D.D.DMgMn[1] ≈ -0.35582434974845456
        @test domain1D.D.DFeFe[1] ≈ 0.24039414025060646
        @test domain1D.D.DFeMg[1] ≈ -0.06679903115794805
        @test domain1D.D.DFeMn[1] ≈ -2.437175247864287
        @test domain1D.D.DMnMn[1] ≈ 3.015178710791253
        @test domain1D.D.DMnMg[1] ≈ -0.001981077957486304
        @test domain1D.D.DMnFe[1] ≈ -0.00523218102384065

        @test domainSph.D.DMgMg[1] ≈ domain1D.D.DMgMg[1]
        @test domainSph.D.DMgFe[1] ≈ domain1D.D.DMgFe[1]
        @test domainSph.D.DMgMn[1] ≈ domain1D.D.DMgMn[1]
        @test domainSph.D.DFeFe[1] ≈ domain1D.D.DFeFe[1]
        @test domainSph.D.DFeMg[1] ≈ domain1D.D.DFeMg[1]
        @test domainSph.D.DFeMn[1] ≈ domain1D.D.DFeMn[1]
        @test domainSph.D.DMnMn[1] ≈ domain1D.D.DMnMn[1]
        @test domainSph.D.DMnMg[1] ≈ domain1D.D.DMnMg[1]
        @test domainSph.D.DMnFe[1] ≈ domain1D.D.DMnFe[1]
    end

    @testset "1D diffusion" begin

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

        IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)

        T = 900u"°C"
        P = 0.6u"GPa"

        domain1D = Domain(IC1D, T, P)

        sol = simulate(domain1D; progress=false, abstol=1e-6,reltol=1e-6, save_everystep=false, save_start=false)

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

        ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lr, tfinal)

        T = 900u"°C"
        P = 0.6u"GPa"

        domainSph = Domain(ICSph, T, P)

        sol = simulate(domainSph; progress=false, abstol=1e-6,reltol=1e-6, save_everystep=false, save_start=false)

        @test norm(sum.(sol.u[end][:,1] .+ sol.u[end][:,2] .+ sol.u[end][:,3])) ≈ 20.272049461615836 rtol=1e-5
    end

    @testset "2D Diffusion" begin

        root = @__DIR__
        path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
        path_2D_Fe = joinpath(root, "Data", "2D", "Xalm_LR.txt")
        path_2D_Mn = joinpath(root, "Data", "2D", "Xsps_LR.txt")
        path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

        CMg = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
        CFe = DelimitedFiles.readdlm(path_2D_Fe, '\t', '\n', header=false)
        CMn = DelimitedFiles.readdlm(path_2D_Mn, '\t', '\n', header=false)
        grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

        Lx = 900.0u"µm"
        Ly = 900.0u"µm"
        tfinal = 1.0u"Myr"
        T = 900u"°C"
        P = 0.6u"GPa"

        IC2D = InitialConditions2D(CMg, CFe, CMn, Lx, Ly, tfinal; grt_boundary = grt_boundary)
        domain2D = Domain(IC2D, T, P)

        sol = simulate(domain2D; save_everystep=false, save_start=false)

        @test norm(sol.u[end][:,:,1]) ≈ 12.783354423875512 rtol=1e-5
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

        IC3D = InitialConditions3D(Mg0, Fe0, Mn0, Lx, Ly, Lz, tfinal; grt_boundary = grt_boundary)
        domain3D = Domain(IC3D, T, P)

        sol = simulate(domain3D; save_everystep=false, save_start=false);

        @test norm(sol.u[end][:,:,:,1]) ≈ 371.1775756471261
    end


    @testset "Callback update D0" begin

        root = @__DIR__
        path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

        data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

        Mg0 = reverse(data[1:size(data,1)÷2, 4])
        Fe0 = reverse(data[1:size(data,1)÷2, 2])
        Mn0 = reverse(data[1:size(data,1)÷2, 3])
        Ca0 = reverse(data[1:size(data,1)÷2, 5])
        distance = data[1:size(data,1)÷2, 1]
        Lx = Lr = (data[end,1] - data[1,1])u"µm"
        tfinal = 1u"Myr"
        T = [850u"°C", 600u"°C"]
        P = [0.5u"GPa", 0.3u"GPa"]

        ICSph = InitialConditionsSpherical(Mg0, Fe0, Mn0, Lr, tfinal)
        IC1D = InitialConditions1D(Mg0, Fe0, Mn0, Lx, tfinal)

        time_update = [0u"Myr", 1u"Myr"]

        domainSph = Domain(ICSph, T, P, time_update)
        domain1D = Domain(IC1D, T, P, time_update)

        @test domainSph.D0[1] ≈ 153548.37186274922
        @test domain1D.D0[1] ≈ 153548.37186274922

        path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
        path_2D_Fe = joinpath(root, "Data", "2D", "Xalm_LR.txt")
        path_2D_Mn = joinpath(root, "Data", "2D", "Xsps_LR.txt")
        path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

        CMg = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
        CFe = DelimitedFiles.readdlm(path_2D_Fe, '\t', '\n', header=false)
        CMn = DelimitedFiles.readdlm(path_2D_Mn, '\t', '\n', header=false)
        grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

        Lx = 900.0u"µm"
        Ly = 900.0u"µm"

        IC2D = InitialConditions2D(CMg, CFe, CMn, Lx, Ly, tfinal; grt_boundary = grt_boundary)
        domain2D = Domain(IC2D, T, P, time_update)

        @test domain2D.D0[1] ≈ 153548.37186274922

        @unpack time_update_ad = domainSph

        update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

        sol_sph = simulate(domainSph; callback=update_diffusion_coef_call, progress=false, abstol=1e-6,reltol=1e-6)

        update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

        sol_1D = simulate(domain1D; callback=update_diffusion_coef_call, progress=false, abstol=1e-6,reltol=1e-6)

        @unpack time_update_ad = domain2D
        update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef)

        sol_2D = simulate(domain2D; callback=update_diffusion_coef_call,save_everystep=false, save_start=false)

        T=600  # in °C
        P=3  # in kbar
        D0 = zeros(Float64, 4)
        DiffusionGarnet.D_ini!(D0,T,P)

        @test sol_1D.prob.p.domain.D0[1] ≈ D0[1]
        @test sol_sph.prob.p.domain.D0[1] ≈ D0[1]
        @test sol_2D.prob.p.domain.D0[1] ≈ D0[1]
    end

    @testset "Callback output" begin

        root = @__DIR__
        path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

        data_1D = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

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

        path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
        path_2D_Fe = joinpath(root, "Data", "2D", "Xalm_LR.txt")
        path_2D_Mn = joinpath(root, "Data", "2D", "Xsps_LR.txt")

        Mg_2D = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
        Fe_2D = DelimitedFiles.readdlm(path_2D_Fe, '\t', '\n', header=false)
        Mn_2D = DelimitedFiles.readdlm(path_2D_Mn, '\t', '\n', header=false)
        Lx = 900.0u"µm"
        Ly = 900.0u"µm"

        IC2D = InitialConditions2D(Mg_2D, Fe_2D, Mn_2D, Lx, Ly, tfinal)
        domain2D = Domain(IC2D, T, P)

        time_save = [0u"Myr", 2u"Myr", 5u"Myr", 15u"Myr"]

        save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain1D.t_charact, save_data)

        sol_1D = simulate(domain1D; callback=save_data_callback, path_save="Grt_1D.h5", abstol=1e-6,reltol=1e-6)

        save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain1D.t_charact, save_data)

        sol_sph = simulate(domainSph; callback=save_data_callback, path_save="Grt_Sph.h5", abstol=1e-6,reltol=1e-6)

        save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain2D.t_charact, save_data)

        sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Grt_2D.h5", progress=false, save_everystep=false)

        h5open("Grt_1D.h5", "r") do file
            @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == convert(Array{Float32}, IC1D.CMg0)
            @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == convert(Array{Float32}, IC1D.CFe0)
            @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == convert(Array{Float32}, IC1D.CMn0)
            @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == convert(Array{Float32}, replace!((1 .- IC1D.CMg0 .- IC1D.CFe0 .- IC1D.CMn0), 1=>0))
            @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == convert(Array{Float32}, sol_1D.u[end][:,1])
            @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == convert(Array{Float32}, sol_1D.u[end][:,2])
            @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == convert(Array{Float32}, sol_1D.u[end][:,3])
        end

        h5open("Grt_Sph.h5", "r") do file
            @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == convert(Array{Float32}, ICSph.CMg0)
            @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == convert(Array{Float32}, ICSph.CFe0)
            @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == convert(Array{Float32}, ICSph.CMn0)
            @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == convert(Array{Float32}, replace!((1 .- ICSph.CMg0 .- ICSph.CFe0 .- ICSph.CMn0), 1=>0))
            @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == convert(Array{Float32}, sol_sph.u[end][:,1])
            @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == convert(Array{Float32}, sol_sph.u[end][:,2])
            @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == convert(Array{Float32}, sol_sph.u[end][:,3])
        end

        h5open("Grt_2D.h5", "r") do file
            @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == convert(Array{Float32}, IC2D.CMg0)
            @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == convert(Array{Float32}, IC2D.CFe0)
            @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == convert(Array{Float32}, IC2D.CMn0)
            @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == convert(Array{Float32}, replace!((1 .- IC2D.CMg0 .- IC2D.CFe0 .- IC2D.CMn0), 1=>0))
            @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == convert(Array{Float32}, sol_2D.u[end][:,:,1])
            @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == convert(Array{Float32}, sol_2D.u[end][:,:,2])
            @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == convert(Array{Float32}, sol_2D.u[end][:,:,3])
        end

        # delete files "Grt_1D.h5" and "Grt_Sph.h5" if it exists
        if isfile("Grt_1D.h5")
            rm("Grt_1D.h5")
        end

        if isfile("Grt_Sph.h5")
            rm("Grt_Sph.h5")
        end

        if isfile("Grt_2D.h5")
            rm("Grt_2D.h5")
        end

        save_data_callback = PresetTimeCallback(ustrip.(time_save) ./ domain2D.t_charact, save_data_paraview)

        sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Grt_2D.h5", progress=false, save_everystep=false, save_start=false)

        h5open("Grt_2D.h5", "r") do file
            @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"])' == convert(Array{Float32}, IC2D.CMg0)
            @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"])' == convert(Array{Float32}, IC2D.CFe0)
            @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"])' == convert(Array{Float32}, IC2D.CMn0)
            @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"])' == convert(Array{Float32},replace!((1 .- IC2D.CMg0 .- IC2D.CFe0 .- IC2D.CMn0), 1=>0))
            @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"])' == convert(Array{Float32}, sol_2D.u[end][:,:,1])
            @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"])' == convert(Array{Float32}, sol_2D.u[end][:,:,2])
            @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"])' == convert(Array{Float32}, sol_2D.u[end][:,:,3])
        end

        # delete files if it exists
        if isfile("Grt_2D.h5")
            rm("Grt_2D.h5")
        end

        if isfile("Grt_2D.xdmf")
            rm("Grt_2D.xdmf")
        end

    end
end
