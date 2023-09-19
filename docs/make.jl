# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter

makedocs(;
         sitename = "DiffusionGarnet.jl",
         modules  = [DiffusionGarnet],
         pages=[
                "Home" => "index.md"
                "Reference" => "reference.md"
               ],
)

deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl", devbranch = "main"
)