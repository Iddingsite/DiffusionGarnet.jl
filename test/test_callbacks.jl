
@testset "Callback update D0" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    CMg0 = reverse(data[1:size(data,1)÷2, 4])
    CFe0 = reverse(data[1:size(data,1)÷2, 2])
    CMn0 = reverse(data[1:size(data,1)÷2, 3])
    CCa0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lx = Lr = (data[end,1] - data[1,1])u"µm"
    tfinal = 1u"Myr"
    T = [850u"°C", 600u"°C"]
    P = [0.5u"GPa", 0.3u"GPa"]

    ICSph = ICSphMajor(;CMg0, CFe0, CMn0, Lr, tfinal)
    IC1D = IC1DMajor(;CMg0, CFe0, CMn0, Lx, tfinal)

    time_update = [0u"Myr", 1u"Myr"]

    domainSph = Domain(ICSph, T, P, time_update)
    domain1D = Domain(IC1D, T, P, time_update)

    @test domainSph.D0[1] ≈ 153548.37186274922
    @test domain1D.D0[1] ≈ 153548.37186274922

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

    IC2D = IC2DMajor(;CMg0, CFe0, CMn0, Lx, Ly, tfinal, grt_boundary)
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

    DiffusionGarnet.D_update!(D0,T_K,P_kbar,D0_data)

    @test sol_1D.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_sph.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_2D.prob.p.domain.D0[1] ≈ D0[1]
end

@testset "Callback output" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data_1D = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    CMg0 = data_1D[:, 4]
    CFe0 = data_1D[:, 2]
    CMn0 = data_1D[:, 3]
    CCa0 = data_1D[:, 5]
    distance = data_1D[:, 1]
    Lx = (data_1D[end,1] - data_1D[1,1])u"µm"
    tfinal = 15u"Myr"

    IC1D = IC1DMajor(;CMg0, CFe0, CMn0, Lx, tfinal)
    ICSph = ICSphMajor(;CMg0, CFe0, CMn0, Lr=Lx, tfinal)

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

    IC2D = IC2DMajor(;CMg0 = Mg_2D, CFe0 = Fe_2D, CMn0 = Mn_2D, Lx, Ly, tfinal)
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

