
@inline function column_to_row(array)
    permutedims(array, reverse(1:ndims(array)))
end

function hdf5_initial_conditions(IC::InitialConditions1D, Domain::Domain1D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["LengthX(µm)"] = IC.Lx
        attributes(g)["Dx(µm)"] = IC.Δx
        attributes(g)["Nx"] = IC.nx
        attributes(g)["TotalTime(Myr)"] = IC.tfinal
        attributes(g)["Coordinates"] = "1D Cartesian"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        create_group(t0, "Mg")
        create_group(t0, "Fe")
        create_group(t0, "Mn")
        create_group(t0, "Ca")

        # describe type of data
        attributes(t0["Mg"])["DataType"] = "Scalar"
        attributes(t0["Fe"])["DataType"] = "Scalar"
        attributes(t0["Mn"])["DataType"] = "Scalar"
        attributes(t0["Ca"])["DataType"] = "Scalar"

        attributes(t0["Mg"])["Center"] = "Node"
        attributes(t0["Fe"])["Center"] = "Node"
        attributes(t0["Mn"])["Center"] = "Node"
        attributes(t0["Ca"])["Center"] = "Node"

        t0["Mg"]["Mg"] = IC.CMg0
        t0["Fe"]["Fe"] = IC.CFe0
        t0["Mn"]["Mn"] = IC.CMn0
        t0["Ca"]["Ca"] = 1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0
    end
end

function hdf5_initial_conditions(IC::InitialConditionsSpherical, Domain::DomainSpherical, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["Radius(µm)"] = IC.Lr
        attributes(g)["Dr(µm)"] = IC.Δr
        attributes(g)["Nr"] = IC.nr
        attributes(g)["TotalTime(Myr)"] = IC.tfinal
        attributes(g)["Coordinates"] = "Spherical"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        create_group(t0, "Mg")
        create_group(t0, "Fe")
        create_group(t0, "Mn")
        create_group(t0, "Ca")

        # describe type of data
        attributes(t0["Mg"])["DataType"] = "Scalar"
        attributes(t0["Fe"])["DataType"] = "Scalar"
        attributes(t0["Mn"])["DataType"] = "Scalar"
        attributes(t0["Ca"])["DataType"] = "Scalar"

        attributes(t0["Mg"])["Center"] = "Node"
        attributes(t0["Fe"])["Center"] = "Node"
        attributes(t0["Mn"])["Center"] = "Node"
        attributes(t0["Ca"])["Center"] = "Node"

        t0["Mg"]["Mg"] = IC.CMg0
        t0["Fe"]["Fe"] = IC.CFe0
        t0["Mn"]["Mn"] = IC.CMn0
        t0["Ca"]["Ca"] = 1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0
    end
end


function hdf5_initial_conditions(IC::InitialConditions2D, Domain::Domain2D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["LengthX(µm)"] = IC.Lx
        attributes(g)["LengthY(µm)"] = IC.Ly
        attributes(g)["Dx(µm)"] = IC.Δx
        attributes(g)["Dy(µm)"] = IC.Δy
        attributes(g)["Nx"] = IC.nx
        attributes(g)["Ny"] = IC.ny
        attributes(g)["TotalTime(Myr)"] = IC.tfinal
        attributes(g)["Coordinates"] = "2D Cartesian"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        attributes(t0)["Temperature(°C)"] = Domain.T[1]
        attributes(t0)["Pressure(GPa)"] = Domain.P[1]
        create_group(t0, "Mg")
        create_group(t0, "Fe")
        create_group(t0, "Mn")
        create_group(t0, "Ca")

        # describe type of data
        attributes(t0["Mg"])["DataType"] = "Scalar"
        attributes(t0["Fe"])["DataType"] = "Scalar"
        attributes(t0["Mn"])["DataType"] = "Scalar"
        attributes(t0["Ca"])["DataType"] = "Scalar"

        attributes(t0["Mg"])["Center"] = "Node"
        attributes(t0["Fe"])["Center"] = "Node"
        attributes(t0["Mn"])["Center"] = "Node"
        attributes(t0["Ca"])["Center"] = "Node"

        t0["Mg"]["Mg"] = IC.CMg0
        t0["Fe"]["Fe"] = IC.CFe0
        t0["Mn"]["Mn"] = IC.CMn0
        t0["Ca"]["Ca"] = 1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0
    end
end

function hdf5_initial_conditions_paraview(IC::InitialConditions2D, Domain::Domain2D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["LengthX(µm)"] = IC.Lx
        attributes(g)["LengthY(µm)"] = IC.Ly
        attributes(g)["Dx(µm)"] = IC.Δx
        attributes(g)["Dy(µm)"] = IC.Δy
        attributes(g)["Nx"] = IC.nx
        attributes(g)["Ny"] = IC.ny
        attributes(g)["TotalTime(Myr)"] = IC.tfinal
        attributes(g)["Coordinates"] = "2D Cartesian"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        attributes(t0)["Temperature(°C)"] = Domain.T[1]
        attributes(t0)["Pressure(GPa)"] = Domain.P[1]
        create_group(t0, "Mg")
        create_group(t0, "Fe")
        create_group(t0, "Mn")
        create_group(t0, "Ca")

        # describe type of data
        attributes(t0["Mg"])["DataType"] = "Scalar"
        attributes(t0["Fe"])["DataType"] = "Scalar"
        attributes(t0["Mn"])["DataType"] = "Scalar"
        attributes(t0["Ca"])["DataType"] = "Scalar"

        attributes(t0["Mg"])["Center"] = "Node"
        attributes(t0["Fe"])["Center"] = "Node"
        attributes(t0["Mn"])["Center"] = "Node"
        attributes(t0["Ca"])["Center"] = "Node"

        t0["Mg"]["Mg"] = column_to_row(IC.CMg0)
        t0["Fe"]["Fe"] = column_to_row(IC.CFe0)
        t0["Mn"]["Mn"] = column_to_row(IC.CMn0)
        t0["Ca"]["Ca"] = column_to_row(1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0)
    end
end

function view_u(u::T1) where {T1 <: AbstractArray{<:Real, 2}}

    CMg = @view u[:,1]
    CFe = @view u[:,2]
    CMn = @view u[:,3]

    return CMg, CFe, CMn
end


function view_u(u::T1) where {T1 <: AbstractArray{<:Real, 3}}

    CMg = @view u[:,:,1]
    CFe = @view u[:,:,2]
    CMn = @view u[:,:,3]

    return CMg, CFe, CMn
end

function view_u(u::T1) where {T1 <: AbstractArray{<:Real, 4}}

    CMg = @view u[:,:,:,1]
    CFe = @view u[:,:,:,2]
    CMn = @view u[:,:,:,3]

    return CMg, CFe, CMn
end


function hdf5_timestep(u, dt, tcurrent, path_hdf5)

    CMg, CFe, CMn = view_u(u)

    h5open(path_hdf5, "r+") do file

        # output the number of group in the HDF5
        n = length((file["Diffusion_Grt"]))

        t = create_group(file, "Diffusion_Grt/t$(lpad(string(n), 4, "0"))") # create a group
        attributes(t)["Time(Myr)"] = tcurrent
        attributes(t)["CurrentDt(Myr)"] = dt
        create_group(t, "Mg")
        create_group(t, "Fe")
        create_group(t, "Mn")
        create_group(t, "Ca")

        # describe type of data
        attributes(t["Mg"])["DataType"] = "Scalar"
        attributes(t["Fe"])["DataType"] = "Scalar"
        attributes(t["Mn"])["DataType"] = "Scalar"
        attributes(t["Ca"])["DataType"] = "Scalar"

        attributes(t["Mg"])["Center"] = "Node"
        attributes(t["Fe"])["Center"] = "Node"
        attributes(t["Mn"])["Center"] = "Node"
        attributes(t["Ca"])["Center"] = "Node"

        t["Mg"]["Mg"] = CMg
        t["Fe"]["Fe"] = CFe
        t["Mn"]["Mn"] = CMn
        t["Ca"]["Ca"] = replace!((1 .- CMg .- CFe .- CMn), 1=>0)
    end
end

function hdf5_timestep_paraview(u, dt, tcurrent, path_hdf5)

    CMg, CFe, CMn = view_u(u)

    h5open(path_hdf5, "r+") do file

        # output the number of group in the HDF5
        n = length((file["Diffusion_Grt"]))

        t = create_group(file, "Diffusion_Grt/t$(lpad(string(n), 4, "0"))") # create a group
        attributes(t)["Time(Myr)"] = tcurrent
        attributes(t)["CurrentDt(Myr)"] = dt
        create_group(t, "Mg")
        create_group(t, "Fe")
        create_group(t, "Mn")
        create_group(t, "Ca")

        # describe type of data
        attributes(t["Mg"])["DataType"] = "Scalar"
        attributes(t["Fe"])["DataType"] = "Scalar"
        attributes(t["Mn"])["DataType"] = "Scalar"
        attributes(t["Ca"])["DataType"] = "Scalar"

        attributes(t["Mg"])["Center"] = "Node"
        attributes(t["Fe"])["Center"] = "Node"
        attributes(t["Mn"])["Center"] = "Node"
        attributes(t["Ca"])["Center"] = "Node"

        t["Mg"]["Mg"] = column_to_row(CMg)
        t["Fe"]["Fe"] = column_to_row(CFe)
        t["Mn"]["Mn"] = column_to_row(CMn)
        t["Ca"]["Ca"] = column_to_row(replace!((1 .- CMg .- CFe .- CMn), 1=>0))
    end
end


function save_data(integrator)

    @unpack IC, t_charact = integrator.p.domain
    @unpack path_save = integrator.p

    if integrator.t ≠ 0.0
        hdf5_timestep(integrator.u, integrator.dt * t_charact, integrator.t * t_charact, path_save)
        println("Data saved at $(round((integrator.t * t_charact), digits=2)) Myr.")
    elseif integrator.t == 0.0
        hdf5_initial_conditions(IC, integrator.p.domain, path_save)
        println("Data saved at $(round((integrator.t * t_charact), digits=2)) Myr.")
    end

end