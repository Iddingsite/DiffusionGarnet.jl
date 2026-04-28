
@testset "Callback update D0 Major" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    CMg0 = reverse(data[1:size(data,1)÷2, 4])
    CFe0 = reverse(data[1:size(data,1)÷2, 2])
    CMn0 = reverse(data[1:size(data,1)÷2, 3])
    CCa0 = reverse(data[1:size(data,1)÷2, 5])
    distance = data[1:size(data,1)÷2, 1]
    Lx = Lr = (data[end,1] - data[1,1])u"µm"
    tfinal = 0.001u"Myr"
    T = [850u"°C", 600u"°C"]
    P = [0.5u"GPa", 0.3u"GPa"]

    ICSph = ICSphMajor(;CMg0, CFe0, CMn0, Lr, tfinal)
    IC1D = IC1DMajor(;CMg0, CFe0, CMn0, Lx, tfinal)

    time_update = [0u"Myr", tfinal]

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

    #now 3D
    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")
    file = jldopen(path_3D, "r")
    Mg0 = file["Mg0"]
    Fe0 = file["Fe0"]
    Mn0 = file["Mn0"]
    close(file)

    Lz = 900.0u"µm"

    IC3D = IC3DMajor(;CMg0 = Mg0, CFe0 = Fe0, CMn0 = Mn0, Lx, Ly, Lz, tfinal)
    domain3D = Domain(IC3D, T, P, time_update)

    @test domain3D.D0[1] ≈ 153548.37186274922

    (; time_update_ad) = domainSph

    update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef, save_positions=(false,false))

    sol_sph = simulate(domainSph; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false, abstol=1e-6,reltol=1e-6)

    sol_1D = simulate(domain1D; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false, abstol=1e-6,reltol=1e-6)

    sol_2D = simulate(domain2D; callback=update_diffusion_coef_call,save_everystep=false, save_start=false, progress=false)

    sol_3D = simulate(domain3D; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false)

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

    D0[1] = ustrip(u"µm^2/Myr", compute_D(Grt_Mg, T = T_K, P = P_kbar))
    D0[2] = ustrip(u"µm^2/Myr", compute_D(Grt_Fe, T = T_K, P = P_kbar))
    D0[3] = ustrip(u"µm^2/Myr", compute_D(Grt_Mn, T = T_K, P = P_kbar))
    D0[4] = 0.5 * D0[2]  # DCa = 0.5 * DFe

    @test sol_1D.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_1D.prob.p.domain.D0[2] ≈ D0[2]
    @test sol_1D.prob.p.domain.D0[3] ≈ D0[3]
    @test sol_1D.prob.p.domain.D0[4] ≈ D0[4]
    @test sol_sph.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_sph.prob.p.domain.D0[2] ≈ D0[2]
    @test sol_sph.prob.p.domain.D0[3] ≈ D0[3]
    @test sol_sph.prob.p.domain.D0[4] ≈ D0[4]
    @test sol_2D.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_2D.prob.p.domain.D0[2] ≈ D0[2]
    @test sol_2D.prob.p.domain.D0[3] ≈ D0[3]
    @test sol_2D.prob.p.domain.D0[4] ≈ D0[4]
    @test sol_3D.prob.p.domain.D0[1] ≈ D0[1]
    @test sol_3D.prob.p.domain.D0[2] ≈ D0[2]
    @test sol_3D.prob.p.domain.D0[3] ≈ D0[3]
    @test sol_3D.prob.p.domain.D0[4] ≈ D0[4]
end

@testset "Callback update D Trace" begin

    nx = 400
    L = 0.8u"mm"
    tfinal = 0.01u"Myr"
    time_update = [0u"Myr", tfinal]
    Grt_REE = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_REE)

    # Gaussian initial condition
    C = zeros(nx)
    C[1:nx÷2] .= 4.5
    C[nx÷2+1:end] .= 1.2

    T = [850u"°C", 600u"°C"]
    P = [0.5u"GPa", 0.3u"GPa"]

    # 1D case
    IC1D = IC1DTrace(;C0=C, D=D, Lx=L, tfinal=0.1u"Myr")
    domain1D = Domain(IC1D, T, P, time_update)
    r = range(0u"mm", length=nx, stop=L)
    ICsph = ICSphTrace(;C0=C, D=D, Lr=L, r=r, tfinal=0.1u"Myr")
    domainsph = Domain(ICsph, T, P, time_update)


    root = @__DIR__
    path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
    path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

    CMg0 = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

    # multiply by 1000 to simulate ppm
    CMg0 .= CMg0 .* 1e5

    Lx = 10.0u"µm"
    Ly = 10.0u"µm"


    IC2D = IC2DTrace(;C0=CMg0, D=D, Lx, Ly, tfinal, grt_boundary)
    domain2D = Domain(IC2D, T, P, time_update)

    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")

    file = jldopen(path_3D, "r")
    Mg0 = file["Mg0"]
    grt_boundary = file["grt_boundary"]
    close(file)


    Lx = 9000.0u"µm"
    Ly = 9000.0u"µm"
    Lz = 9000.0u"µm"

    Mg0 .= Mg0 .* 1e5

    IC3D = IC3DTrace(;C0=Mg0, D=D, Lx, Ly, Lz, tfinal, grt_boundary)
    domain3D = Domain(IC3D, T, P, time_update)

    (; time_update_ad) = domainsph
    update_diffusion_coef_call = PresetTimeCallback(time_update_ad, update_diffusion_coef, save_positions=(false,false))

    sol_1D = simulate(domain1D; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false, abstol=1e-6,reltol=1e-6)
    sol_sph = simulate(domainsph; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false, abstol=1e-6,reltol=1e-6)
    sol_2D = simulate(domain2D; callback=update_diffusion_coef_call,save_everystep=false, save_start=false, progress=false)
    sol_3D = simulate(domain3D; callback=update_diffusion_coef_call, save_everystep=false, save_start=false, progress=false)

    # compute expected D at T=600°C and P=3kbar

    D_val = ustrip(u"µm^2/Myr", compute_D(D, T = T[end], P = P[end]))


    @test sol_1D.prob.p.domain.D[1] ≈ D_val
    @test sol_sph.prob.p.domain.D[1] ≈ D_val
    @test sol_2D.prob.p.domain.D[1] ≈ D_val
    @test sol_3D.prob.p.domain.D[1] ≈ D_val
end


@testset "Callback output Major" begin

    root = @__DIR__
    path_1D = joinpath(root, "Data", "1D", "Data_Grt_1D.txt")

    data_1D = DelimitedFiles.readdlm(path_1D, '\t', '\n', header=true)[1]

    CMg0 = data_1D[:, 4]
    CFe0 = data_1D[:, 2]
    CMn0 = data_1D[:, 3]
    CCa0 = data_1D[:, 5]
    distance = data_1D[:, 1]
    Lx = (data_1D[end,1] - data_1D[1,1])u"µm"
    tfinal = 0.01u"Myr"

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

    #now 3D
    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")
    file = jldopen(path_3D, "r")
    Mg0 = file["Mg0"]
    Fe0 = file["Fe0"]
    Mn0 = file["Mn0"]
    close(file)

    Lz = 900.0u"µm"

    IC3D = IC3DMajor(;CMg0 = Mg0, CFe0 = Fe0, CMn0 = Mn0, Lx, Ly, Lz, tfinal)
    domain3D = Domain(IC3D, T, P)

    time_save = [0u"Myr", 0.001u"Myr", 0.005u"Myr", tfinal]

    save_data_callback = PresetTimeCallback(ustrip.(time_save), save_data, save_positions=(false,false))

    sol_1D = simulate(domain1D; callback=save_data_callback, path_save="Grt_1D.h5", abstol=1e-6,reltol=1e-6)

    sol_sph = simulate(domainSph; callback=save_data_callback, path_save="Grt_Sph.h5", abstol=1e-6,reltol=1e-6)

    sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Grt_2D.h5", progress=false, save_everystep=false)

    sol_3D = simulate(domain3D; callback=save_data_callback, path_save="Grt_3D.h5", progress=false, save_everystep=false)


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

    # 3D
    h5open("Grt_3D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"]) == convert(Array{Float32}, IC3D.CMg0)
        @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"]) == convert(Array{Float32}, IC3D.CFe0)
        @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"]) == convert(Array{Float32}, IC3D.CMn0)
        @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"]) == convert(Array{Float32}, replace!((1 .- IC3D.CMg0 .- IC3D.CFe0 .- IC3D.CMn0), 1=>0))
        @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"]) == convert(Array{Float32}, sol_3D.u[end][:,:,:,1])
        @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"]) == convert(Array{Float32}, sol_3D.u[end][:,:,:,2])
        @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"]) == convert(Array{Float32}, sol_3D.u[end][:,:,:,3])
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

    save_data_callback = PresetTimeCallback(ustrip.(time_save), save_data_paraview, save_positions=(false,false))

    sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Grt_2D.h5", progress=false, save_everystep=false, save_start=false)

    sol_3D = simulate(domain3D; callback=save_data_callback, path_save="Grt_3D.h5", progress=false, save_everystep=false, save_start=false)

    h5open("Grt_2D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"])' ≈ convert(Array{Float32}, IC2D.CMg0)
        @test read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"])' ≈ convert(Array{Float32}, IC2D.CFe0)
        @test read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"])' ≈ convert(Array{Float32}, IC2D.CMn0)
        @test read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"])' ≈ convert(Array{Float32},replace!((1 .- IC2D.CMg0 .- IC2D.CFe0 .- IC2D.CMn0), 1=>0))
        @test read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"])' ≈ convert(Array{Float32}, sol_2D.u[end][:,:,1])
        @test read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"])' ≈ convert(Array{Float32}, sol_2D.u[end][:,:,2])
        @test read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"])' ≈ convert(Array{Float32}, sol_2D.u[end][:,:,3])
    end

    h5open("Grt_3D.h5", "r") do file
        array = read(file["Diffusion_Grt"]["t0000"]["Mg"]["Mg"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, IC3D.CMg0)
        array = read(file["Diffusion_Grt"]["t0000"]["Fe"]["Fe"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, IC3D.CFe0)
        array = read(file["Diffusion_Grt"]["t0000"]["Mn"]["Mn"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, IC3D.CMn0)
        array = read(file["Diffusion_Grt"]["t0000"]["Ca"]["Ca"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, replace!((1 .- IC3D.CMg0 .- IC3D.CFe0 .- IC3D.CMn0), 1=>0))
        array = read(file["Diffusion_Grt"]["t0003"]["Mg"]["Mg"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, sol_3D.u[end][:,:,:,1])
        array = read(file["Diffusion_Grt"]["t0003"]["Fe"]["Fe"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, sol_3D.u[end][:,:,:,2])
        array = read(file["Diffusion_Grt"]["t0003"]["Mn"]["Mn"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, sol_3D.u[end][:,:,:,3])
    end

    # delete files if it exists
    if isfile("Grt_2D.h5")
        rm("Grt_2D.h5")
    end

    if isfile("Grt_2D.xdmf")
        rm("Grt_2D.xdmf")
    end

    if isfile("Grt_3D.h5")
        rm("Grt_3D.h5")
    end

    if isfile("Grt_3D.xdmf")
        rm("Grt_3D.xdmf")
    end

end

@testset "Callback output Trace" begin
    nx = 400
    L = 0.8u"mm"
    Grt_REE = Garnet.Grt_REE_Bloch2020_fast
    D = SetChemicalDiffusion(Grt_REE)

    # Gaussian initial condition
    C = zeros(nx)
    C[1:nx÷2] .= 4.5
    C[nx÷2+1:end] .= 1.2

    # 1D case
    IC1D = IC1DTrace(;C0=C, D=D, Lx=L, tfinal=0.1u"Myr")
    domain1D = Domain(IC1D, 800u"°C", 1u"GPa")
    r = range(0u"mm", length=nx, stop=L)
    ICsph = ICSphTrace(;C0=C, D=D, Lr=L, r=r, tfinal=0.1u"Myr")
    domainsph = Domain(ICsph, 800u"°C", 1u"GPa")


    root = @__DIR__
    path_2D_Mg = joinpath(root, "Data", "2D", "Xprp_LR.txt")
    path_2D_grt = joinpath(root, "Data", "2D", "Contour_LR.txt")

    CMg0 = DelimitedFiles.readdlm(path_2D_Mg, '\t', '\n', header=false)
    grt_boundary = DelimitedFiles.readdlm(path_2D_grt, '\t', '\n', header=false)

    # multiply by 1000 to simulate ppm
    CMg0 .= CMg0 .* 1e5

    Lx = 10.0u"µm"
    Ly = 10.0u"µm"
    tfinal = 0.1u"Myr"

    IC2D = IC2DTrace(;C0=CMg0, D=D, Lx, Ly, tfinal, grt_boundary)
    domain2D = Domain(IC2D, 800u"°C", 1u"GPa")

    path_3D = joinpath(root, "Data", "3D", "3D_data.jld2")

    file = jldopen(path_3D, "r")
    Mg0 = file["Mg0"]
    grt_boundary = file["grt_boundary"]
    close(file)

    Lx = 9000.0u"µm"
    Ly = 9000.0u"µm"
    Lz = 9000.0u"µm"
    T = 800u"°C"
    P = 1u"GPa"

    Mg0 .= Mg0 .* 1e5

    IC3D = IC3DTrace(;C0=Mg0, D=D, Lx, Ly, Lz, tfinal, grt_boundary)
    domain3D = Domain(IC3D, T, P)

    time_save = [0u"Myr", 0.01u"Myr", tfinal]

    save_data_callback = PresetTimeCallback(ustrip.(time_save), save_data, save_positions=(false,false))

    sol_1D = simulate(domain1D; callback=save_data_callback, path_save="Trace_1D.h5", abstol=1e-6,reltol=1e-6, save_everystep=false)
    sol_sph = simulate(domainsph; callback=save_data_callback, path_save="Trace_Sph.h5", abstol=1e-6,reltol=1e-6, save_everystep=false)

    time_save = [0u"Myr", 0.005u"Myr", 0.01u"Myr"]

    save_data_callback = PresetTimeCallback(ustrip.(time_save), save_data, save_positions=(false,false))

    sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Trace_2D.h5", progress=false, save_everystep=false)
    sol_3D = simulate(domain3D; callback=save_data_callback, path_save="Trace_3D.h5", progress=false, save_everystep=false)

    h5open("Trace_1D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["C"]["C"]) ≈ convert(Array{Float32}, IC1D.C0)
        @test read(file["Diffusion_Grt"]["t0002"]["C"]["C"]) ≈ convert(Array{Float32}, sol_1D.u[end])
    end

    h5open("Trace_Sph.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["C"]["C"]) ≈ convert(Array{Float32}, ICsph.C0)
        @test read(file["Diffusion_Grt"]["t0002"]["C"]["C"]) ≈ convert(Array{Float32}, sol_sph.u[end])
    end

    h5open("Trace_2D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["C"]["C"]) ≈ convert(Array{Float32}, IC2D.C0)
        @test read(file["Diffusion_Grt"]["t0002"]["C"]["C"]) ≈ convert(Array{Float32}, sol_2D.u[end])
    end

    h5open("Trace_3D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["C"]["C"]) ≈ convert(Array{Float32}, IC3D.C0)
        @test read(file["Diffusion_Grt"]["t0002"]["C"]["C"]) ≈ convert(Array{Float32}, sol_3D.u[end])
    end

    # delete files if it exists
    if isfile("Trace_1D.h5")
        rm("Trace_1D.h5")
    end

    if isfile("Trace_Sph.h5")
        rm("Trace_Sph.h5")
    end

    if isfile("Trace_2D.h5")
        rm("Trace_2D.h5")
    end

    if isfile("Trace_3D.h5")
        rm("Trace_3D.h5")
    end

    save_data_callback = PresetTimeCallback(ustrip.(time_save), save_data_paraview, save_positions=(false,false))

    sol_2D = simulate(domain2D; callback=save_data_callback, path_save="Trace_2D.h5", progress=false, save_everystep=false, save_start=false)

    h5open("Trace_2D.h5", "r") do file
        @test read(file["Diffusion_Grt"]["t0000"]["C"]["C"])' ≈ convert(Array{Float32}, IC2D.C0)
        @test read(file["Diffusion_Grt"]["t0002"]["C"]["C"])' ≈ convert(Array{Float32}, sol_2D.u[end])
    end

    # do the same for 3D
    sol_3D = simulate(domain3D; callback=save_data_callback, path_save="Trace_3D.h5", progress=false, save_everystep=false, save_start=false)

    h5open("Trace_3D.h5", "r") do file
        array = read(file["Diffusion_Grt"]["t0000"]["C"]["C"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, IC3D.C0)
        array = read(file["Diffusion_Grt"]["t0002"]["C"]["C"])
        @test permutedims(array,1:ndims(array)) ≈ convert(Array{Float32}, sol_3D.u[end])
    end

    # delete files if it exists
    if isfile("Trace_2D.h5")
        rm("Trace_2D.h5")
    end

    if isfile("Trace_2D.xdmf")
        rm("Trace_2D.xdmf")
    end

    if isfile("Trace_3D.h5")
        rm("Trace_3D.h5")
    end

    if isfile("Trace_3D.xdmf")
        rm("Trace_3D.xdmf")
    end
end
