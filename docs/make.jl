# Inside make.jl
push!(LOAD_PATH,"../src/")
using DiffusionGarnet
using Documenter

makedocs(;
         sitename = "DiffusionGarnet.jl",
         modules  = [DiffusionGarnet],
         pages=[
                "Home" => "index.md"
               ],
)

deploydocs(;
    repo="github.com/Iddingsite/DiffusionGarnet.jl.git", devbranch = "main"
)