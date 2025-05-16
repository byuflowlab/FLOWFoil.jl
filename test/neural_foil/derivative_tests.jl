@testset "Check NeuralFoil Derivatives" begin
    function wrapfun(var)
        x, y = at.naca4(var[1], var[2], var[3])
        coordinates = [x y]
        flow_angles = range(-5, 5, 3)
        reynolds = 1e6
        mach = 0.0
        model_size = "xlarge"

        outputs_ff = analyze(
            coordinates,
            flow_angles;
            method=NeuralFoil(reynolds, mach; model_size=model_size),
        )

        return [outputs_ff.cl; outputs_ff.cd; outputs_ff.cm; outputs_ff.confidence]
    end

    nacaparam = [2.0, 4.0, 12.0]
    adjac = ForwardDiff.jacobian(wrapfun, nacaparam)
    fdjac = FiniteDiff.finite_difference_jacobian(wrapfun, nacaparam)

    @test isapprox(adjac, fdjac, atol=1e-7)
end
