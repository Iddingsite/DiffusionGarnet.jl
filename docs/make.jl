# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter

makedocs(;
         sitename = "DiffusionGarnet.jl",
         pages=[
                "Home" => "index.md"
                "Tutorials" => Any["Diffusion in 1D Cartesian coordinates" => "diffusion_1D.md",
                "Diffusion in spherical coordinates" => "diffusion_spherical.md",
                "Diffusion in 2D Cartesian coordinates on CPU" => "diffusion_2D.md"
                ]
                "Callbacks" => Any["Updating PT conditions" => "updating_PT.md",
                "Saving output" => "saving_output.md"]
                "Reference" => "reference.md"
               ],
)

deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl", devbranch = "main"
)