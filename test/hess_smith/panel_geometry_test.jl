@testset "Hess Smith" begin

    x = [1.0, 0.5, 0.0, 0.5, 1.0]
    y = [0.0, -0.5, 0.0, 0.5, 0.0]
    xy = [x y]

    # Generate Panel Geometry
    panel_geometry = FLOWFoil.generate_panel_geometry(FLOWFoil.HessSmith(), xy)
    @test panel_geometry.npanels == 4
    @test panel_geometry.panel_center == [0.75 -0.25; 0.25 -0.25; 0.25 0.25; 0.75 0.25]
    @test isapprox(panel_geometry.panel_length, ones(4) * 0.70710678, atol=1e-4)
    @test isapprox(panel_geometry.sine_vector, [-0.7071067811865475, 0.7071067811865475, 0.7071067811865475, -0.7071067811865475], atol=1e-6)
    @test isapprox(panel_geometry.cosine_vector, [-0.7071067811865475, -0.7071067811865475, 0.7071067811865475, 0.7071067811865475], atol=1e-6)
end