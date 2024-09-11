@testset "Axisymmetric System Assembly Tests" begin

    #---------------------------------#
    #           AXISYMMETRIC          #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01] .+ 1.0
    coordinates = [x y]
    pt = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false])

    # Generate Panel Geometry
    panels = generate_panels(pt, coordinates)

    # Generate Influence Mesh
    mesh = generate_mesh(pt, panels)

    # Assemble Linear System
    system = generate_inviscid_system(pt, panels, mesh)

    # check that it's not invertable
    @test !isapprox(LinearAlgebra.det(system.A), 0.0)

    # # check that the kutta condition is in the right place
    # TODO: this needs to be updated to a new test for the updated kutta condition implementation
    # @test system.A[end, 1] == 1.0
    # @test system.A[end, end - 1] == 1.0
    # @test system.A[end, end] == 0.0
    # @test all(system.A[end, 2:(end - 2)] .== 0.0)

    # TODO: need to test boundary condition matrix

    # - Multi Body Test - #

    x1 = [0.0; 0.5; 1.0]
    r1 = [0.0; 0.5; 0.0]

    x2 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r2 = [1.5; 1.0; 1.5; 2.0; 1.5]

    coordinates = ([x1 r1], [x2 r2])
    pt = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [true, false])
    panel_array = generate_panels(pt, coordinates)
    mesh = generate_mesh(pt, panel_array)
    system = generate_inviscid_system(pt, panel_array, mesh)

    # @test system.A[:, end] == [0; 0; 1; 1; 1; 1; 0.0]
    # @test system.A[end, :] == [0, 0, 1, 0, 0, 1, 0.0]
    @test size(system.A) == (5, 5)
    @test isapprox(
        system.b,
        [
            -sqrt(2) / 2
            -sqrt(2) / 2
            sqrt(2)
            sqrt(2) / 2
            -sqrt(2) / 2
        ],
    )

    #TODO add back diagonal correction test.

end


