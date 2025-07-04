@testset "Periodic System Assembly Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#

    # 4 panel diamond test
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [0.0; -1.0; 0.0; 1.0; 0.0]
    coordinates = [x y]

    #test planar (method 1) and cascade(method 2)
    system_geometry_method_1 = FLOWFoil.Martensen(false, 0.0, 0.0, 30.0, false)
    system_geometry_method_2 = FLOWFoil.Martensen(true, 2.0, 0.0, 30.0, false)

    panels_1 = FLOWFoil.generate_panel_geometry(system_geometry_method_1, coordinates)
    panels_2 = FLOWFoil.generate_panel_geometry(system_geometry_method_2, coordinates)

    mesh_1 = FLOWFoil.generate_system_geometry(system_geometry_method_1, panels_1)
    mesh_2 = FLOWFoil.generate_system_geometry(system_geometry_method_2, panels_2)

    # Assemble Linear System
    system_1_matrices = FLOWFoil.generate_system_matrices(
        system_geometry_method_1, panels_2, mesh_1
    )
    system_2_matrices = FLOWFoil.generate_system_matrices(
        system_geometry_method_2, panels_2, mesh_2
    )

    #test A matrices
    @test isapprox(
        system_1_matrices.A,
        [
            -1.108732318685386 0.1909859317102744 -0.19098593171027442
            0.1909859317102744 -0.5 0.054366159342693045
            -0.1909859317102744 0.054366159342692996 -0.5
        ],
        atol=0.001,
    )
    @test isapprox(
        system_2_matrices.A,
        [
            6.000055797672284 0.0 0.0
            0.0 -0.5 -3.500027898836142
            0.0 -3.500027898836142 -0.5
        ],
        atol=0.001,
    )

    #test right hand side vectors
    @test isapprox(
        system_1_matrices.b,
        [0.894423 1.18672e-6; 0.447215 -0.89442658; -0.44721 -0.89442776],
        atol=0.001,
    )
    @test isapprox(
        system_2_matrices.b,
        [0.894423 1.18672e-6; 0.447215 -0.89442658; -0.44721 -0.89442776],
        atol=0.001,
    )
end

