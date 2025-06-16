@testset "NeuralFoil Tests" begin

    # set operating conditions
    flow_angles = -5:1:15
    reynolds = 1e6
    model_size = "xlarge"

    x, y = FLOWFoil.AirfoilTools.naca4()

    nf_outputs = nf.get_aero_from_coordinates(
        reverse([x y]; dims=1), flow_angles, reynolds; model_size=model_size
    )

    outputs = FLOWFoil.analyze(
        [x y], flow_angles; method=NeuralFoil(reynolds; model_size=model_size)
    )

    @test isapprox(outputs.cl, nf_outputs.cl[:, 1])
    @test isapprox(outputs.cd, nf_outputs.cd[:, 1])
    @test isapprox(outputs.cm, nf_outputs.cm[:, 1])
    @test isapprox(outputs.analysis_confidence, nf_outputs.analysis_confidence[:, 1], atol=1e-5)
end
