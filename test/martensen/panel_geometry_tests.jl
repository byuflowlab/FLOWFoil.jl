@testset "Periodic Panel Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#
    # 4 panel diamond test
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    z = [0.0; -1.0; 0.0; 1.0; 0.0]
    coordinates = [x z]
    panel_geometry_method = FLOWFoil.Martensen(false, 0.0, 0.0, 30.0, 100.0, false)
    panels = FLOWFoil.generate_panel_geometry(
        panel_geometry_method,
        coordinates
    )
    @test panels.npanels == 4 #test number of panels
    @test isapprox(panels.panel_angle, [-2.03444, -4.24874, 1.10715, -1.10715], atol = 0.001) #test panel angle
    @test isapprox(panels.delta_angle, [-1.01722, 1.570795, 1.570795, -0.553575], atol = 0.001) #test delta angle
end
