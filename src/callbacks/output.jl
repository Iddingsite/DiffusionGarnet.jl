function hdf5_initial_conditions(IC::InitialConditions1D, Domain::Domain1D, path_hdf5::String)
    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        # Set attributes for group g
        attrs_values = (
            ("LengthX(µm)", IC.Lx),
            ("Dx(µm)", IC.Δx),
            ("Nx", IC.nx),
            ("TotalTime(Myr)", IC.tfinal),
            ("Coordinates", "1D Cartesian"),
            ("CharacteristicLength", Domain.L_charact),
            ("CharacteristicDiffusionCoefficient", Domain.D_charact),
            ("CharacteristicTime", Domain.t_charact)
        )
        for (attr, value) in attrs_values
            attributes(g)[attr] = value
        end

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca")
        for group in groups
            grp = create_group(t0, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t0["Mg"]["Mg"] = convert(Array{Float32}, IC.CMg0)
        t0["Fe"]["Fe"] = convert(Array{Float32}, IC.CFe0)
        t0["Mn"]["Mn"] = convert(Array{Float32}, IC.CMn0)
        t0["Ca"]["Ca"] = convert(Array{Float32}, replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
    end
end


function hdf5_initial_conditions(IC::InitialConditionsSpherical, Domain::DomainSpherical, path_hdf5::String)
    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        # Set attributes for group g
        attrs_values = (
            ("Radius(µm)", IC.Lr),
            ("Dr(µm)", IC.Δr),
            ("Nr", IC.nr),
            ("TotalTime(Myr)", IC.tfinal),
            ("Coordinates", "Spherical"),
            ("CharacteristicLength", Domain.L_charact),
            ("CharacteristicDiffusionCoefficient", Domain.D_charact),
            ("CharacteristicTime", Domain.t_charact)
        )
        for (attr, value) in attrs_values
            attributes(g)[attr] = value
        end

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca")
        for group in groups
            grp = create_group(t0, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t0["Mg"]["Mg"] = convert(Array{Float32}, IC.CMg0)
        t0["Fe"]["Fe"] = convert(Array{Float32}, IC.CFe0)
        t0["Mn"]["Mn"] = convert(Array{Float32}, IC.CMn0)
        t0["Ca"]["Ca"] = convert(Array{Float32}, replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
    end
end


function hdf5_initial_conditions(IC::InitialConditions2D, Domain::Domain2D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        # Set attributes for group g
        attrs_values = (
            ("LengthX(µm)", IC.Lx),
            ("LengthY(µm)", IC.Ly),
            ("Dx(µm)", IC.Δx),
            ("Dy(µm)", IC.Δy),
            ("Nx", IC.nx),
            ("Ny", IC.ny),
            ("TotalTime(Myr)", IC.tfinal),
            ("Coordinates", "2D Cartesian"),
            ("CharacteristicLength", Domain.L_charact),
            ("CharacteristicDiffusionCoefficient", Domain.D_charact),
            ("CharacteristicTime", Domain.t_charact)
        )
        for (attr, value) in attrs_values
            attributes(g)[attr] = value
        end

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        attributes(t0)["Temperature(°C)"] = Domain.T[1]
        attributes(t0)["Pressure(GPa)"] = Domain.P[1]

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca", "GrtPosition", "GrtBoundary")
        for group in groups
            grp = create_group(t0, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t0["Mg"]["Mg"] = convert(Array{Float32}, IC.CMg0)
        t0["Fe"]["Fe"] = convert(Array{Float32}, IC.CFe0)
        t0["Mn"]["Mn"] = convert(Array{Float32}, IC.CMn0)
        t0["Ca"]["Ca"] = convert(Array{Float32}, replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
        t0["GrtPosition"]["GrtPosition"] = convert(Array{Int32},IC.grt_position)
        t0["GrtBoundary"]["GrtBoundary"] = convert(Array{Int32},IC.grt_boundary)
    end
end

function hdf5_initial_conditions(IC::InitialConditions3D, Domain::Domain3D, path_hdf5)

    h5open(path_hdf5, "w") do file
        g = create_group(file, "Diffusion_Grt") # create a group

        # Set attributes for group g
        attrs_values = (
            ("LengthX(µm)", IC.Lx),
            ("LengthY(µm)", IC.Ly),
            ("LengthZ(µm)", IC.Lz),
            ("Dx(µm)", IC.Δx),
            ("Dy(µm)", IC.Δy),
            ("Dz(µm)", IC.Δz),
            ("Nx", IC.nx),
            ("Ny", IC.ny),
            ("Nz", IC.nz),
            ("TotalTime(Myr)", IC.tfinal),
            ("Coordinates", "3D Cartesian"),
            ("CharacteristicLength", Domain.L_charact),
            ("CharacteristicDiffusionCoefficient", Domain.D_charact),
            ("CharacteristicTime", Domain.t_charact)
        )
        for (attr, value) in attrs_values
            attributes(g)[attr] = value
        end

        t0 = create_group(file, "Diffusion_Grt/t$(lpad("0", 4, "0"))") # create a group
        attributes(t0)["Time(Myr)"] = 0
        attributes(t0)["Temperature(°C)"] = Domain.T[1]
        attributes(t0)["Pressure(GPa)"] = Domain.P[1]

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca", "GrtPosition", "GrtBoundary")
        for group in groups
            grp = create_group(t0, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t0["Mg"]["Mg"] = convert(Array{Float32}, IC.CMg0)
        t0["Fe"]["Fe"] = convert(Array{Float32}, IC.CFe0)
        t0["Mn"]["Mn"] = convert(Array{Float32}, IC.CMn0)
        t0["Ca"]["Ca"] = convert(Array{Float32}, replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
        t0["GrtPosition"]["GrtPosition"] = convert(Array{Int32},IC.grt_position)
        t0["GrtBoundary"]["GrtBoundary"] = convert(Array{Int32},IC.grt_boundary)
    end
end


function view_u(u::T1) where {T1 <: AbstractArray{<:Real, N}} where N

    indices = ntuple(d -> :, N-1)

    CMg = @view u[indices..., 1]
    CFe = @view u[indices..., 2]
    CMn = @view u[indices..., 3]

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

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca")
        for group in groups
            grp = create_group(t, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t["Mg"]["Mg"] = convert(Array{Float32}, CMg)
        t["Fe"]["Fe"] = convert(Array{Float32}, CFe)
        t["Mn"]["Mn"] = convert(Array{Float32}, CMn)
        t["Ca"]["Ca"] = convert(Array{Float32}, replace!((1 .- CMg .- CFe .- CMn), 1=>0))
    end
end

"""
    save_data(integrator)

Callback function used to save major element compositions to an HDF5 file at a specific timestep.
"""
function save_data(integrator)

    @unpack IC, t_charact = integrator.p.domain
    @unpack path_save = integrator.p

    if integrator.t ≠ 0.0
        hdf5_timestep(integrator.u, integrator.dt * t_charact, integrator.t * t_charact, path_save)
        @info "Data saved at $(round((integrator.t * t_charact), digits=2)) Myr."
    elseif integrator.t == 0.0
        hdf5_initial_conditions(IC, integrator.p.domain, path_save)
        @info "Data saved at $(round((integrator.t * t_charact), digits=2)) Myr."
    end
end