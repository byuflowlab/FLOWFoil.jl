@testset "Periodic Mesh Tests" begin

    # - Very Basic Test - #

    x = [0.0; 1.0; 2.0]
    z = [0.0; 1.0; 2.0]
    coordinates = [x z]
    panels = generate_panels(
        PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0), coordinates
    )
    mesh = generate_mesh(PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0), panels)

    @test mesh.nbodies == 1
    @test mesh.panel_indices == [1:2]
    @test mesh.x == [0.0 -1.0; 1.0 0.0]
    @test mesh.y == [0.0 -1.0; 1.0 0.0]
end
