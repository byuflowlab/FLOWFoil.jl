@testset "Post Processing Tests" begin

    # - Very Basic Test - #

    x, z = AirfoilTools.naca4()
    coordinates = [x z]
    outputs = analyze(coordinates, [0.0; 5.0]; method=Mfoil())

    set_coordinates(reverse(x), reverse(z))
    cl0, cm0 = solve_alpha(0.0)
    cl5, cm5 = solve_alpha(5.0)

    @test isapprox(outputs.cl[1], cl0, atol=1e-5)
    @test isapprox(outputs.cm[1], cm0, atol=1e-3)
    @test isapprox(outputs.cl[2], cl5, atol=1e-5)
    @test isapprox(outputs.cm[2], cm5, atol=1e-3)
end
