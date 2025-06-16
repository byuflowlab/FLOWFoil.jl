@testset "LegacyXfoil Tests" begin

    # set operating conditions
    flow_angles = -5:1:15
    reynolds = 1e6

    x, y = FLOWFoil.AirfoilTools.naca4()

    clxf, cdxf, _, cmxf, convxf = xf.alpha_sweep(
        reverse(x), reverse(y), collect(flow_angles), reynolds
    )

    outputs = FLOWFoil.analyze([x y], flow_angles; method=LegacyXfoil(reynolds))

    @test isapprox(outputs.cl, clxf)
    @test isapprox(outputs.cd, cdxf)
    @test isapprox(outputs.cm, cmxf)
    @test outputs.converged == convxf
end
