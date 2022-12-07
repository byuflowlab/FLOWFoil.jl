@testset "Planar Mesh Tests" begin

    # - Very Basic Test - #

    x = [1.0; 0.0; 1.0]
    z = [0.0; 0.0; eps()]
    coordinates = [x z]
    panels = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    mesh, TEmesh = generate_mesh(PlanarProblem(Vortex(Linear()), Dirichlet()), panels)

    @test mesh.panel_indices == [1:2]
    @test mesh.node_indices == [1:3]
    @test mesh.chord == 1.0
    @test mesh.panel_length == [1.0; 1.0]
    @test mesh.r1 == [0.0 1.0; 1.0 0.0; eps() 1.0]
    @test mesh.lnr1 == zeros(3, 2)
    @test mesh.r1normal == [0.0 -eps(); 0.0 0.0; -eps() 0.0]
    @test mesh.r1tangent == [0.0 1.0; 1.0 0.0; 0.0 1.0]
    @test mesh.theta1 == [pi 0.0; 0.0 pi; pi 0.0]
    @test mesh.r2 == [1.0 eps(); 0.0 1.0; 1.0 0.0]
    @test mesh.lnr2 == zeros(3, 2)
    @test mesh.theta2 == [pi 0.0; 0.0 pi; pi 0.0]

    @test TEmesh.blunt_te == Bool[0]
    @test TEmesh.trailing_edge_gap == [0.0]
    @test TEmesh.tdp == [1.0]
    @test TEmesh.txp == [0.0]
end

@testset "Axisymmetric Mesh Tests" begin

    # - Very Basic Test - #

    x = [0.0; 1.0; 2.0]
    z = [0.0; 1.0; 2.0]
    coordinates = [x z]
    panels = generate_panels(
        AxisymmetricProblem([true], Vortex(Constant()), Neumann()), coordinates
    )
    mesh = generate_mesh(AxisymmetricProblem([true], Vortex(Constant()), Neumann()), panels)

    @test mesh.nbodies == 1
    @test mesh.panel_indices == [1:2]
    @test mesh.x == [0.0 (0.5 - 1.5)/1.5; 1.0/0.5 0.0]
    @test mesh.r == [1.0 0.5/1.5; 1.5/0.5 1.0]
    @test mesh.m == [1.0 0.6; 0.6 1.0]
end
