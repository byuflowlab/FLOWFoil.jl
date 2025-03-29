@testset "Post Processing Tests" begin

    x = [1.0, 0.5, 0.0, 0.5, 1.0]
    y = [0.0, -0.5, 0.0, 0.5, 0.0]
    xy = [x y]
    flow_angles = [4.0]

    # Generate Panel Geometry
    panel_geometry = FLOWFoil.generate_panel_geometry(FLOWFoil.HessSmith(), xy)
    system_geometry = FLOWFoil.generate_system_geometry(FLOWFoil.HessSmith(), panel_geometry)
    system_matrices = FLOWFoil.generate_system_matrices(FLOWFoil.HessSmith(), panel_geometry, system_geometry)
    strengths = FLOWFoil.solve(FLOWFoil.HessSmith(), system_matrices)
    outputs = FLOWFoil.post_process(FLOWFoil.HessSmith(), panel_geometry, system_geometry, strengths, flow_angles)

    @test isapprox(outputs.tangential_velocities, [[-1.2255216557014643 -1.0282209987511828 1.7933098485248005 1.596009191574519]], atol=1e-4)
    @test isapprox(outputs.surface_pressures, [[-0.5019033285932584 -0.057238422272879896 -2.215960212816043 -1.5472453395903494]], atol=1e-4)

end
