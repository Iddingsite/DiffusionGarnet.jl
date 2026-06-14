using DiffusionGarnet
using Test
using Unitful

# Only meaningful when the active backend is Float32 + 3D.
# On other backends the IC/Domain Float32 tests in test_discretisation_major.jl
# already verify struct-level type propagation.
@testset "Type agnosticism (Float32 end-to-end)" begin
    if occursin("Float32", DiffusionGarnet.backend) && occursin("3D", DiffusionGarnet.backend)
        n = 7
        CMg0 = fill(0.2f0, n, n, n)
        CFe0 = fill(0.6f0, n, n, n)
        CMn0 = fill(0.1f0, n, n, n)

        IC = IC3DMajor(; CMg0, CFe0, CMn0,
                       Lx = 1000.0u"µm", Ly = 1000.0u"µm", Lz = 1000.0u"µm",
                       tfinal = 0.001u"Myr")
        domain = Domain(IC, 900u"°C", 0.6u"GPa")

        @test eltype(domain.u0) == Float32
        @test domain.Δxad_ isa Float32

        sol = simulate(domain; save_everystep = false, save_start = false, progress = false)

        @test eltype(sol.u[end]) == Float32
    else
        @test_skip "backend $(DiffusionGarnet.backend) is not Float32+3D — skipping end-to-end Float32 test"
    end
end
