@testset "Problem Definition Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # Very Basic Test
    coordinates = [[1], [1]]
    p = define_problem(
        PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates, 0.0, -1.0, -1.0
    )
    @test p.nbody == 2
    @test p.angle_of_attack == [0.0]
    @test p.reynolds == [-1.0]
    @test p.mach == [-1.0]
    @test p.viscous == false
    @test p.method == PlanarProblem(Vortex(Linear()), Dirichlet())
end
