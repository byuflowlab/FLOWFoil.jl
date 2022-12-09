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
    @test p.flow_angle == [0.0]
    @test p.reynolds == [-1.0]
    @test p.mach == [-1.0]
    @test p.viscous == false
    @test p.method == PlanarProblem(Vortex(Linear()), Dirichlet())

    #---------------------------------#
    #          Axisymmetric           #
    #---------------------------------#
    # Very Basic Test
    coordinates = [[1], [1]]
    p = define_problem(
        AxisymmetricProblem(Vortex(Constant()), Neumann(), [false, true]),
        coordinates,
        0.0,
        -1.0,
        -1.0,
    )
    @test p.nbody == 2
    @test p.flow_angle == [0.0]
    @test p.reynolds == [-1.0]
    @test p.mach == [-1.0]
    @test p.viscous == false
    @test all(p.method.body_of_revolution .== [false, true])

    #---------------------------------#
    #            Periodic             #
    #---------------------------------#
    # Very Basic Test
    coordinates = [[1], [1]]
    p = define_problem(
        PeriodicProblem(Vortex(Constant()), Neumann(), [1.0]),
        coordinates,
        0.0,
        -1.0,
        -1.0,
    )
    @test p.nbody == 2
    @test p.flow_angle == [0.0]
    @test p.reynolds == [-1.0]
    @test p.mach == [-1.0]
    @test p.viscous == false
end
