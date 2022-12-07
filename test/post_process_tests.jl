@testset "Post Processing Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    x, z = x, z = FLOWFoil.naca4()
    coordinates = [x z]
    polar = solve(
        coordinates, [0.0; 5.0]; method=PlanarProblem(Vortex(Linear()), Dirichlet())
    )

    Xfoil.set_coordinates(reverse(x), reverse(z))
    cl0, cm0 = Xfoil.solve_alpha(0.0)
    cl5, cm5 = Xfoil.solve_alpha(5.0)

    @test isapprox(polar.lift[1], cl0, atol=1e-5)
    @test isapprox(polar.moment[1], cm0, atol=1e-5)
    @test isapprox(polar.lift[2], cl5, atol=1e-5)
    @test isapprox(polar.moment[2], cm5, atol=1e-5)
end
