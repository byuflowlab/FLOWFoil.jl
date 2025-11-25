@testset "Check Lewis Method Derivatives" begin
    function wrapfun(var)
        x, y = AirfoilTools.naca4(var[1], var[2], var[3])

        angles_of_attack=range(0.0,5.0)
        outputs = analyze(x, y.+1.0, angles_of_attack, method=Lewis(body_of_revolution=[false]))

        return [outputs.cl; outputs.cd]
    end

    nacaparam = [2.0, 4.0, 12.0]
    adjac = ForwardDiff.jacobian(wrapfun, nacaparam)
    fdjac = FiniteDifferences.jacobian(central_fdm(5, 1), wrapfun, nacaparam)[1]

    println(maximum(adjac.-fdjac))
    @test isapprox(adjac, fdjac, atol=1e-7)
end
