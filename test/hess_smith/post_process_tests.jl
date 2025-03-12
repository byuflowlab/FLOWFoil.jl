
@testset "Post Processing Tests" begin

    x = [1.0, 0.5, 0.0, 0.5, 1.0]
    y = [0.0, -0.5, 0.0, 0.5, 0.0]
    xy = [x y]

    # Generate Panel Geometry
    panel_geometry = FLOWFoil.generate_panel_geometry(FLOWFoil.HessSmith(), xy)
    system_geometry = FLOWFoil.generate_system_geometry(FLOWFoil.HessSmith(), panel_geometry)
    system_matrices = FLOWFoil.generate_system_matrices(FLOWFoil.HessSmith(), panel_geometry, system_geometry)
    strengths = FLOWFoil.solve(FLOWFoil.HessSmith(), system_matrices)
    outputs = FLOWFoil.post_process(FLOWFoil.HessSmith(), panel_geometry, system_geometry, strengths, flow_angles)

    @test 1 == 1 

end
