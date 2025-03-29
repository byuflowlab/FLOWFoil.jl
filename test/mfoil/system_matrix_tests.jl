@testset "Planar System Assembly Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01]
    coordinates = [x y]
    pt = Mfoil()

    # Generate Panel Geometry
    panels = FLOWFoil.generate_panel_geometry(pt, coordinates)

    # Generate Influence Mesh
    mesh = FLOWFoil.generate_system_geometry(pt, panels)

    # Assemble Linear System
    system = FLOWFoil.generate_system_matrices(pt, panels, mesh)

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
    panels = FLOWFoil.generate_panel_geometry(Mfoil(), coordinates)
    mesh = FLOWFoil.generate_system_geometry(Mfoil(), panels)
    system = FLOWFoil.generate_system_matrices(pt, panels, mesh)

    @test all(
        system.b .==
        [0.5 0.0; 0.0 0.0; -0.5 0.0; 0.0 1.0; 0.0 1.5; 0.0 2.0; 0.0 0.0; 0.0 0.0],
    )
end

