@testset "Panel Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # Very Basic Test
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panels = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    @test panels.npanels == 3
    @test panels.panel_edges[1, :, :] == [1.0 0.0; 2.0 0.0]
    @test panels.panel_edges[2, :, :] == [2.0 0.0; 3.0 0.0]
    @test panels.panel_edges[3, :, :] == [3.0 0.0; 4.0 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_vector == [1.0 0.0; 1.0 0.0; 1.0 0.0]

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

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#
    # Very Basic Test
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panels = generate_panels(
        PeriodicProblem(Vortex(Constant()), Neumann(), [1.0]), coordinates
    )
    @test panels.npanels == 3
    @test panels.panel_center == [1.5 0.0; 2.5 0.0; 3.5 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_normal == [0.0 1.0; 0.0 1.0; 0.0 1.0]
    @test panels.panel_angle == zeros(3)
end
