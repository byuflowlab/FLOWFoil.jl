@testset "Planar Panel Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#

    # - Very Basic Test - #
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panel_array = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    panels = panel_array[1]
    @test panels.npanels == 3
    @test panels.panel_edges[1, :, :] == [1.0 0.0; 2.0 0.0]
    @test panels.panel_edges[2, :, :] == [2.0 0.0; 3.0 0.0]
    @test panels.panel_edges[3, :, :] == [3.0 0.0; 4.0 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_vector == [1.0 0.0; 1.0 0.0; 1.0 0.0]

    # - Multi-body Tests - #
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    z1 = [0.0; -0.5; 0.0; 0.5; 0.0]

    x2 = 1.25 .+ [1.0; 0.5; 0.0; 0.5; 1.0]
    z2 = [0.0; -0.5; 0.0; 0.5; 0.0]

    coordinates = ([x1 z1], [x2 z2])
    panel_array = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)

    @test panel_array[1].npanels == 4
    @test panel_array[1].panel_edges[1, :, :] == [1.0 0.0; 0.5 -0.5]
    @test panel_array[1].panel_edges[2, :, :] == [0.5 -0.5; 0.0 0.0]
    @test panel_array[1].panel_edges[3, :, :] == [0.0 0.0; 0.5 0.5]
    @test panel_array[1].panel_edges[4, :, :] == [0.5 0.5; 1.0 0.0]
    @test panel_array[1].panel_length == sqrt(0.5) * ones(4)
    @test panel_array[1].panel_vector == [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]

    @test panel_array[2].npanels == 4
    @test panel_array[2].panel_edges[1, :, :] == [1.25 0.0] .+ [1.0 0.0; 0.5 -0.5]
    @test panel_array[2].panel_edges[2, :, :] == [1.25 0.0] .+ [0.5 -0.5; 0.0 0.0]
    @test panel_array[2].panel_edges[3, :, :] == [1.25 0.0] .+ [0.0 0.0; 0.5 0.5]
    @test panel_array[2].panel_edges[4, :, :] == [1.25 0.0] .+ [0.5 0.5; 1.0 0.0]
    @test panel_array[2].panel_length == sqrt(0.5) * ones(4)
    @test panel_array[2].panel_vector == [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]
end

@testset "Axisymmetric Panel Tests" begin
    #---------------------------------#
    #           AXISYMMETRIC          #
    #---------------------------------#
    # Very Basic Test
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panels = generate_panels(
        AxisymmetricProblem(Vortex(Constant()), Neumann(), [true]), coordinates
    )
    @test panels.npanels == 3
    @test panels.panel_center == [1.5 0.0; 2.5 0.0; 3.5 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_normal == [0.0 1.0; 0.0 1.0; 0.0 1.0]
    @test panels.panel_curvature == zeros(3)
    @test panels.panel_angle == zeros(3)
end

@testset "Periodic Panel Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#
    # Very Basic Test
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panels = generate_panels(
        PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0), coordinates
    )
    @test panels.npanels == 3
    @test panels.panel_center == [1.5 0.0; 2.5 0.0; 3.5 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_normal == [0.0 1.0; 0.0 1.0; 0.0 1.0]
    @test panels.panel_angle == zeros(3)
end
