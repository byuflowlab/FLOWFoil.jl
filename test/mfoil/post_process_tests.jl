@testset "Post Processing Tests" begin

    # - Very Basic Test - #

    x, z = x, z = FLOWFoil.naca4()
    coordinates = [x z]
    outputs = analyze(coordinates, [0.0; 5.0]; method=Mfoil())

    Xfoil.set_coordinates(reverse(x), reverse(z))
    cl0, cm0 = Xfoil.solve_alpha(0.0)
    cl5, cm5 = Xfoil.solve_alpha(5.0)

    @test isapprox(outputs.lift[1], cl0, atol=1e-5)
    @test isapprox(outputs.moment[1], cm0, atol=1e-3)
    @test isapprox(outputs.lift[2], cl5, atol=1e-5)
    @test isapprox(outputs.moment[2], cm5, atol=1e-3)
end
