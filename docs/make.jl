# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter

makedocs(;
         sitename = "DiffusionGarnet.jl",
         pages=[
                "Home" => "index.md"
                "Background" => "background.md"
                "Tutorials" => Any["Diffusion in 1D Cartesian coordinates" => "diffusion_1D.md",
                "Diffusion in spherical coordinates" => "diffusion_spherical.md",
                "Diffusion in 2D Cartesian coordinates on CPU" => "diffusion_2D.md",
                "Callbacks" => Any["Updating pressure and temperature conditions" => "updating_PT.md",
                "Saving output as HDF5 files" => "saving_output.md"],
                "Diffuion in 3D Cartesian coordinates on CPU" => "diffusion_3D.md"
                ]
                "List of functions" => "reference.md"
               ],
        format = Documenter.HTML(; mathengine=
        Documenter.MathJax3(Dict(  # use MathJax3 as engine for latex (to be able to reference equations)
            :loader => Dict("load" => ["[tex]/physics"]),
            :tex => Dict(
                "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
                "tags" => "ams",
            )))),
)

deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl", devbranch = "main"
)