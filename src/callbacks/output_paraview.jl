
@inline function column_to_row(array)
    permutedims(array, reverse(1:ndims(array)))
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
        create_group(t0, "GrtPosition")
        create_group(t0, "GrtBoundary")

        # describe type of data
        attributes(t0["Mg"])["DataType"] = "Scalar"
        attributes(t0["Fe"])["DataType"] = "Scalar"
        attributes(t0["Mn"])["DataType"] = "Scalar"
        attributes(t0["Ca"])["DataType"] = "Scalar"
        attributes(t0["GrtPosition"])["DataType"] = "Scalar"
        attributes(t0["GrtBoundary"])["DataType"] = "Scalar"

        attributes(t0["Mg"])["Center"] = "Node"
        attributes(t0["Fe"])["Center"] = "Node"
        attributes(t0["Mn"])["Center"] = "Node"
        attributes(t0["Ca"])["Center"] = "Node"
        attributes(t0["GrtPosition"])["Center"] = "Node"
        attributes(t0["GrtBoundary"])["Center"] = "Node"

        t0["Mg"]["Mg"] = column_to_row(IC.CMg0)
        t0["Fe"]["Fe"] = column_to_row(IC.CFe0)
        t0["Mn"]["Mn"] = column_to_row(IC.CMn0)
        t0["Ca"]["Ca"] = column_to_row(replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
        t0["GrtPosition"]["GrtPosition"] = IC.grt_position
        t0["GrtBoundary"]["GrtBoundary"] = IC.grt_boundary
    end
end

function hdf5_initial_conditions_paraview(IC::InitialConditions3D, Domain::Domain3D, path_hdf5)

  h5open(path_hdf5, "w") do file
      g = create_group(file, "Diffusion_Grt") # create a group

      attributes(g)["LengthX(µm)"] = IC.Lx
      attributes(g)["LengthY(µm)"] = IC.Ly
      attributes(g)["LengthZ(µm)"] = IC.Lz
      attributes(g)["Dx(µm)"] = IC.Δx
      attributes(g)["Dy(µm)"] = IC.Δy
      attributes(g)["Dz(µm)"] = IC.Δz
      attributes(g)["Nx"] = IC.nx
      attributes(g)["Ny"] = IC.ny
      attributes(g)["Nz"] = IC.nz
      attributes(g)["TotalTime(Myr)"] = IC.tfinal
      attributes(g)["Coordinates"] = "3D Cartesian"
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
      create_group(t0, "GrtPosition")
      create_group(t0, "GrtBoundary")

      # describe type of data
      attributes(t0["Mg"])["DataType"] = "Scalar"
      attributes(t0["Fe"])["DataType"] = "Scalar"
      attributes(t0["Mn"])["DataType"] = "Scalar"
      attributes(t0["Ca"])["DataType"] = "Scalar"
      attributes(t0["GrtPosition"])["DataType"] = "Scalar"
      attributes(t0["GrtBoundary"])["DataType"] = "Scalar"

      attributes(t0["Mg"])["Center"] = "Node"
      attributes(t0["Fe"])["Center"] = "Node"
      attributes(t0["Mn"])["Center"] = "Node"
      attributes(t0["Ca"])["Center"] = "Node"
      attributes(t0["GrtPosition"])["Center"] = "Node"
      attributes(t0["GrtBoundary"])["Center"] = "Node"

      t0["Mg"]["Mg"] = column_to_row(IC.CMg0)
      t0["Fe"]["Fe"] = column_to_row(IC.CFe0)
      t0["Mn"]["Mn"] = column_to_row(IC.CMn0)
      t0["Ca"]["Ca"] = column_to_row(replace!((1 .- IC.CMg0 .- IC.CFe0 .- IC.CMn0), 1=>0))
      t0["GrtPosition"]["GrtPosition"] = IC.grt_position
      t0["GrtBoundary"]["GrtBoundary"] = IC.grt_boundary
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

        # Create groups and set attributes for each group
        groups = ("Mg", "Fe", "Mn", "Ca")
        for group in groups
            grp = create_group(t, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        t["Mg"]["Mg"] = column_to_row(CMg)
        t["Fe"]["Fe"] = column_to_row(CFe)
        t["Mn"]["Mn"] = column_to_row(CMn)
        t["Ca"]["Ca"] = column_to_row(replace!((1 .- CMg .- CFe .- CMn), 1=>0))
    end
end


function XMDF_creation(path_hdf5)

    # create xdmf file name by removing everything after the last "." and adding ".xdmf"
    xdmf_file_path = splitext(path_hdf5)[1] * ".xdmf"
    # xdmf file name without path
    xdmf_file_name = split(xdmf_file_path, "/")[end]
    # hdf5 file name without path

    hdf5_file_name = split(path_hdf5, "/")[end]

    hdf = h5open(path_hdf5, "r")
    Grt = hdf["Diffusion_Grt"]

    open(xdmf_file_path, "w") do file
        # write header
        write(file, "<?xml version=\"1.0\" ?>
    <!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>
    <Xdmf Version=\"3.0\">
      <Domain>
        <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">")

        # extract Nz, Ny, Nx
        nz = haskey(attributes(Grt), "Nz") ? read_attribute(Grt, "Nz") : 0
        ny = haskey(attributes(Grt), "Ny") ? read_attribute(Grt, "Ny") : 0
        nx = haskey(attributes(Grt), "Nx") ? read_attribute(Grt, "Nx") : 0

        # extract dz, dy, dx
        dz = haskey(attributes(Grt), "Dz(µm)") ? read_attribute(Grt, "Dz(µm)") : 0
        dy = haskey(attributes(Grt), "Dy(µm)") ? read_attribute(Grt, "Dy(µm)") : 0
        dx = haskey(attributes(Grt), "Dx(µm)") ? read_attribute(Grt, "Dx(µm)") : 0

        # write timesteps
        for(i, t) in enumerate(Grt)
            write(file, "
          <Grid Name=\"mesh$(string(i-1))\" GridType=\"Uniform\">
            <Time Value=\"$(read_attribute(t,"Time(Myr)"))\" />
            <!-- provide Nz, Ny, Nx -->
            <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"$nz $ny $nx\"/>
            <Geometry GeometryType=\"ORIGIN_DXDYDZ\">
              <!-- Oz,Oy,Ox + Dz,Dy,Dx-->
              <DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">
              <!-- where start origins -->
                0.0 0.0 0.0
              </DataItem>
              <DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">
                <!-- dz input.dy input.dx -->
                $dz $dy $dx
              </DataItem>
            </Geometry>")

            timestep = keys(Grt)[i]  # t0000, t0001, t0002, ...

            for j in 1:length(t)
                element = keys(t)[j]  # Ca, Mg, Fe or Mn

                write(file, "
            <Attribute Name=\"$(element)\" AttributeType=\"Scalar\" Center=\"Node\">
              <!-- provide Nz, Ny, Nx -->
              <DataItem Dimensions=\"$nz $ny $nx\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">
                $(hdf5_file_name):/Diffusion_Grt/$(timestep)/$(element)/$(element)
              </DataItem>
            </Attribute>")
            end

            # close timestep
            write(file, "
          </Grid>")
        end

        # write footer of xdmf file
        write(file, "
        </Grid>
      </Domain>
    </Xdmf>")
    end

    close(hdf)
end


function save_data_paraview(integrator)

    @unpack IC, t_charact = integrator.p.domain
    @unpack path_save = integrator.p

    if integrator.t ≠ 0.0
        hdf5_timestep_paraview(integrator.u, integrator.dt * t_charact, integrator.t * t_charact, path_save)
        XMDF_creation(path_save)
        println("Data saved at $(round((integrator.t * t_charact), digits=2)) Myr.")
    elseif integrator.t == 0.0
        hdf5_initial_conditions_paraview(IC, integrator.p.domain, path_save)
        XMDF_creation(path_save)
        println("Data saved at $(round((integrator.t * t_charact), digits=2)) Myr.")
    end

end