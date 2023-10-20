using HDF5

function hdf5_initial_conditions(IC::InitialConditions1D, Domain::Domain1D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["LengthX"] = IC.Lx
        attributes(g)["Dx"] = IC.Δx
        attributes(g)["Nx"] = IC.nx
        attributes(g)["TotalTimeYears"] = IC.tfinal
        attributes(g)["Coordinates"] = "1D Cartesian"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time"] = 0
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

function hdf5_initial_conditions(IC::InitialConditionsSpherical, Domain::DomainSpherical, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["Radius"] = IC.Lr
        attributes(g)["Dr"] = IC.Δr
        attributes(g)["Nr"] = IC.nr
        attributes(g)["TotalTimeYears"] = IC.tfinal
        attributes(g)["Coordinates"] = "Spherical"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time"] = 0
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


function hdf5_initial_conditions(IC::InitialConditions2D, Domain::Domain2D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        attributes(g)["LengthX"] = IC.Lx
        attributes(g)["LengthY"] = IC.Ly
        attributes(g)["Dx"] = IC.Δx
        attributes(g)["Dy"] = IC.Δy
        attributes(g)["Nx"] = IC.nx
        attributes(g)["Nx"] = IC.ny
        attributes(g)["TotalTimeYears"] = IC.tfinal
        attributes(g)["Coordinates"] = "2D Cartesian"
        attributes(g)["CharacteristicLength"] =  Domain.L_charact
        attributes(g)["CharacteristicDiffusionCoefficient"] =  Domain.D_charact
        attributes(g)["CharacteristicTime"] =  Domain.t_charact

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time"] = 0
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

function hdf5_timestep(IC, Domain, dt, tcurrent)

    @unpack path_data = IC

    h5open(path_data, "r+") do file
        t = create_group(file, "Diffusion_Grt/t$(lpad(string(tcurrent), 4, "0"))") # create a group
        attributes(t)["Time"] = tcurrent
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

        t["Mg"]["Mg"] = column_to_row(Domain.CMg)
        t["Fe"]["Fe"] = column_to_row(Domain.CFe)
        t["Mn"]["Mn"] = column_to_row(Domain.CMn)
        t["Ca"]["Ca"] = column_to_row(1 .- Domain.CMg .- Domain.CFe .- Domain.CMn)
    end
end





end