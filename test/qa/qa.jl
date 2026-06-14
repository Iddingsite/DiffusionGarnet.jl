using Aqua
using DiffusionGarnet

@testset "Aqua" begin
    Aqua.test_all(
        DiffusionGarnet;
        ambiguities    = false,          # ParallelStencil macros generate many
        stale_deps     = (ignore = [:CUDA],),   # CUDA loaded via ParallelStencil, not directly
        deps_compat    = (ignore = [:CUDA], check_extras = false),
    )
end
