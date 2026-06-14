using Test

@testset "DiffusionGarnet QA" begin
    include("qa.jl")
    include("type_agnosticism.jl")
    include("jet.jl")
end
