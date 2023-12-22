# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter

makedocs(;
         sitename = "DiffusionGarnet.jl",
         pages=[
                "Home" => "index.md"
                "Tutorials" => Any["Diffusion in 1D" => "diffusion_1D.md"]
                # "Diffusion in 2D" => "diffusion_2D.md",
                # "Diffusion in Sperical Coordinates" => "diffusion_spherical.md"]
                # "Callbacks" => "callbacks.md"
                "Reference" => "reference.md"
               ],
)

deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl", devbranch = "main"
)