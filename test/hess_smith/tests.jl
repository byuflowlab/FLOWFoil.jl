@testset "Hess Smith" begin
    x, y = [diamond]
    xy = [x y]

    method = ff.HessSmith()

    # Generate Panel Geometry
    panel_geometry = ff.generate_panel_geometry(method, xy)

    @test isapprox(panel_geometry.control_points, [])

    # Generate Influence Mesh
    # system_geometry = ff.generate_system_geometry(method, panel_geometry)

    # Assemble Linear System
    # system_matrices = ff.generate_system_matrices(method, panel_geometry, system_geometry)

    # # Solve System
    # strengths = solve(method, system_matrices)

    # # Post Process Solution
    # post_process(method, panel_geometry, system_geometry, strengths, flow_angles)

end
