@testset "Post Processing Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    x, z = x, z = FLOWFoil.naca4()
    coordinates = [x z]
    polar = solve(coordinates, PlanarProblem(Vortex(Linear()), Neumann()))

    Xfoil.set_coordinates(reverse(x), reverse(z))
    cl, cm = Xfoil.solve_alpha(0.0)

    @test isapprox(polar.lift[1], cl, atol=1e-5)
    @test isapprox(polar.moment[1], cm, atol=1e-5)
end
