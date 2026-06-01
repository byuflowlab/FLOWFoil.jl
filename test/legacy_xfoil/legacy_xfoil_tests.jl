@testset "LegacyXfoil Tests" begin

    # set operating conditions
    flow_angles = -5:1:15
    reynolds = 1e6

    x, y = FLOWFoil.AirfoilTools.naca4()

    clxf, cdxf, _, cmxf, convxf = xf.alpha_sweep(
        reverse(x), reverse(y), collect(flow_angles), reynolds
    )

    outputs = FLOWFoil.analyze(
        [x y], flow_angles; method=LegacyXfoil(; reynolds=reynolds)
    )

    @test isapprox(outputs.cl, clxf)
    @test isapprox(outputs.cd, cdxf)
    @test isapprox(outputs.cm, cmxf)
    @test outputs.converged == convxf

    # multi-α does not populate cp / bl
    @test outputs.cp === nothing
    @test outputs.x_cp === nothing
    @test outputs.bl === nothing
end

@testset "LegacyXfoil single-α viscous cp + BL dump" begin
    x, y = FLOWFoil.AirfoilTools.naca4()

    outputs = FLOWFoil.analyze(
        [x y], [5.0]; method=LegacyXfoil(; reynolds=1e6, npan=160)
    )

    @test outputs.cp !== nothing
    @test outputs.x_cp !== nothing
    @test outputs.bl !== nothing
    @test length(outputs.cp) == length(outputs.x_cp)
    @test length(outputs.bl.theta) == length(outputs.bl.s)
    @test all(outputs.bl.theta .>= 0.0)
    @test outputs.converged[1] == true
end

@testset "LegacyXfoil inviscid mode (reynolds=nothing)" begin
    x, y = FLOWFoil.AirfoilTools.naca4()

    outputs = FLOWFoil.analyze([x y], [5.0]; method=LegacyXfoil(; npan=160))

    @test outputs.cd == zeros(1)
    @test outputs.cdp == zeros(1)
    @test outputs.converged == trues(1)
    @test outputs.cp !== nothing
    @test outputs.x_cp !== nothing
    @test outputs.bl === nothing
    @test length(outputs.cp) == length(outputs.x_cp)
end
