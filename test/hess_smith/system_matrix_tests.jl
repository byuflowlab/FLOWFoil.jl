@testset "System Assembly Tests" begin

    # - Very Basic Test - #

    x = [1.0, 0.5, 0.0, 0.5, 1.0]
    y = [0.0, -0.5, 0.0, 0.5, 0.0]
    xy = [x y]

    # Generate Panel Geometry
    panel_geometry = FLOWFoil.generate_panel_geometry(FLOWFoil.HessSmith(), xy)
    system_geometry = FLOWFoil.generate_system_geometry(FLOWFoil.HessSmith(), panel_geometry)
    system_matrices = FLOWFoil.generate_system_matrices(FLOWFoil.HessSmith(), panel_geometry, system_geometry)
    strengths = FLOWFoil.solve(FLOWFoil.HessSmith(), system_matrices)

    # check that the A matrix is not invertable
    @test !isapprox(LinearAlgebra.det(system_matrices.A), 0.0)

    # check that the A matrix assembled correctly
    @test isapprox(system_matrices.A[end, end], 3.946311609806842, atol=1e-6)
    @test isapprox(system_matrices.A[1, 1], 3.1415926535897922, atol=1e-6)

    # check that the B matrix assembled correctly
    @test isapprox(system_matrices.b, 
                [-0.7071067811865475 0.7071067811865475 
                0.7071067811865475 0.7071067811865475 
                0.7071067811865475 -0.7071067811865475 
                -0.7071067811865475 -0.7071067811865475 
                -0.0 1.414213562373095], 
                atol=1e-6
                )

    @test isapprox(strengths,
                [-0.31933685593539957 0.31933685593539957
                0.31933685593539957 0.3193368559353996
                0.3193368559353996 -0.3193368559353996
                -0.31933685593539957 -0.31933685593539957 
                -2.8133181927807724e-17 0.7167267576431025], 
                atol=1e-6
                )
end
