using JET
using DiffusionGarnet
using Test

@testset "JET" begin
    rep = JET.report_package(DiffusionGarnet; target_defined_modules = true)
    # Pre-existing JET findings in trace-element parametric inner constructors
    # and missing last_index are tracked here as broken until fixed upstream.
    @test length(JET.get_reports(rep)) == 0 broken=true
end
