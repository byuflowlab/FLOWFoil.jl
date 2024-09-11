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
