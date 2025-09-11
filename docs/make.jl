# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter
using DocumenterVitepress

makedocs(;
         sitename = "DiffusionGarnet.jl",
         authors="Hugo Dominguez",
         pages=[
                "Home" => "index.md"
                "Background" => "background.md"
                "Tutorials" => Any["Diffusion in 1D Cartesian coordinates" => "diffusion_1D.md",
                "Diffusion in spherical coordinates" => "diffusion_spherical.md",
                "Using different diffusion coefficient datasets" => "comparing_diffcoef.md",
                "Diffusion in 2D Cartesian coordinates on CPU" => "diffusion_2D.md",
                "Callbacks" => Any["Updating pressure and temperature conditions" => "updating_PT.md",
                "Saving output as HDF5 files" => "saving_output.md"],
                "Diffusion in 3D Cartesian coordinates on CPU" => "diffusion_3D_CPU.md",
                "Diffusion in 3D Cartesian coordinates on GPU" => "diffusion_3D_GPU.md"
                ]
                "List of functions" => "reference.md"
               ],
        format = DocumenterVitepress.MarkdownVitepress(; repo="https://github.com/Iddingsite/DiffusionGarnet.jl",
        devbranch = "main",
        devurl = "dev"),
)

DocumenterVitepress.deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl", push_preview=true, target="build", devbranch = "main"
)