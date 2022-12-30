@testset "Planar System Assembly Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01]
    coordinates = [x y]
    pt = PlanarProblem(Vortex(Linear()), Dirichlet())

    # Generate Panel Geometry
    panels = generate_panels(pt, coordinates)

    # Generate Influence Mesh
    mesh, TEmesh = generate_mesh(pt, panels)

    # Assemble Linear System
    system = generate_inviscid_system(pt, panels, mesh, TEmesh)

    # check that it's not invertable
    @test !isapprox(LinearAlgebra.det(system.A), 0.0)

    # check that the kutta condition is in the right place
    @test system.A[end, 1] == 1.0
    @test system.A[end, end - 1] == 1.0
    @test system.A[end, end] == 0.0
    @test all(system.A[end, 2:(end - 2)] .== 0.0)

    # check that the stream function values are in the right places
    @test all(system.A[1:(end - 2), end] .== -1.0)
    @test system.A[end, end] .== 0.0

    # TODO: need to test boundary condition matrix

    # - Multi-body Test - #
    x1 = zeros(3)
    z1 = [-0.5; 0.0; 0.5]

    x2 = [1.0; 1.5; 2.0]
    z2 = zeros(3)

    coordinates = ([x1 z1], [x2 z2])
    panels = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    mesh, TEmesh = generate_mesh(PlanarProblem(Vortex(Linear()), Dirichlet()), panels)
    system = generate_inviscid_system(pt, panels, mesh, TEmesh)

    @test all(
        system.b .==
        [0.5 0.0; 0.0 0.0; -0.5 0.0; 0.0 1.0; 0.0 1.5; 0.0 2.0; 0.0 0.0; 0.0 0.0],
    )
end

@testset "Axisymmetric System Assembly Tests" begin

    #---------------------------------#
    #           AXISYMMETRIC          #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01] .+ 1.0
    coordinates = [x y]
    pt = AxisymmetricProblem(Vortex(Constant()), Neumann(), [false])

    # Generate Panel Geometry
    panels = generate_panels(pt, coordinates)

    # Generate Influence Mesh
    mesh = generate_mesh(pt, panels)

    # Assemble Linear System
    system = generate_inviscid_system(pt, panels, mesh)

    # check that it's not invertable
    @test !isapprox(LinearAlgebra.det(system.A), 0.0)

    # check that the kutta condition is in the right place
    @test system.A[end, 1] == 1.0
    @test system.A[end, end - 1] == 1.0
    @test system.A[end, end] == 0.0
    @test all(system.A[end, 2:(end - 2)] .== 0.0)

    # check that the stream function values are in the right places
    @test all(system.A[1:(end - 1), end] .== 1.0)
    @test all(system.A[end, end] .== 0.0)

    # TODO: need to test boundary condition matrix
    # TODO: multi-body test
end

@testset "Periodic System Assembly Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01] .+ 1.0
    coordinates = [x y]
    pt = PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0)

    # Generate Panel Geometry
    panels = generate_panels(pt, coordinates)

    # Generate Influence Mesh
    mesh = generate_mesh(pt, panels)

    # Assemble Linear System
    system = generate_inviscid_system(pt, panels, mesh)

    # check that it's not invertable
    @test !isapprox(LinearAlgebra.det(system.A), 0.0)

    # check that the kutta condition is in the right place
    @test system.A[end, 1] == 1.0
    @test system.A[end, end - 1] == 1.0
    @test system.A[end, end] == 0.0
    @test all(system.A[end, 2:(end - 2)] .== 0.0)

    # check that the stream function values are in the right places
    @test all(system.A[1:(end - 1), end] .== 1.0)
    @test all(system.A[end, end] .== 0.0)

    # TODO: multi-airfoil test
end
