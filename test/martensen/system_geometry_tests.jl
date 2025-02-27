@testset "Periodic Mesh Tests" begin

    # - Very Basic Test - #

    # 4 panel diamond test
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    z = [0.0; -1.0; 0.0; 1.0; 0.0]
    coordinates = [x z]
    system_geometry_method = FLOWFoil.Martensen(false, 0.0, 0.0, 30.0, 100.0, false)
    panels = FLOWFoil.generate_panel_geometry(
        system_geometry_method,
        coordinates
    )
    mesh = FLOWFOil.generate_system_geometry(
        system_geometry_method,
        panels
    )

    @test mesh.nbodies == 1
    @test mesh.panel_indices == [1:4]
    @test mesh.r_x == [0.0 -0.5 -0.5 0.0; 0.5 0.0 0.0 0.5; 0.5 0.0 0.0 0.5; 0.0 -0.5 -0.5 0.0]
    @test mesh.r_y == [0.0 0.0 1.0 1.0; 0.0 0.0 1.0 1.0; -1.0 -1.0 0.0 0.0; -1.0 -1.0 0.0 0.0]
    @test mesh.r_squared == [0.0 0.25 1.25 1.0; 0.25 0.0 1.0 1.25; 1.25 1.0 0.0 0.25; 1.0 1.25 0.25 0.0]
end
