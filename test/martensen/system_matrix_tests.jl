@testset "Periodic System Assembly Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#

    # 4 panel diamond test
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    z = [0.0; -1.0; 0.0; 1.0; 0.0]
    coordinates = [x z]

    #test planar (method 1) and cascade(method 2)
    system_geometry_method_1 = FLOWFoil.Martensen(false, 0.0, 0.0, 30.0, 100.0)
    system_geometry_method_2 = FLOWFoil.Martensen(true, 1.0, 0.0, 30.0, 100.0)

    panels_1 = FLOWFoil.generate_panel_geometry(
        system_geometry_method_1,
        coordinates
    )
    panels_2 = FLOWFoil.generate_panel_geometry(
        system_geometry_method_2,
        coordinates
    )

    mesh_1 = FLOWFoil.generate_system_geometry(
        system_geometry_method_1,
        panels_1
    )
    mesh_2 = FLOWFoil.generate_system_geometry(
        system_geometry_method_2,
        panels_2
    )

    # Assemble Linear System
    system_1_matrices = FLOWFoil.generate_system_matrices(system_geometry_method_1, panels_2, mesh_1)
    system_2_matrices = FLOWFoil.generate_system_matrices(system_geometry_method_2, panels_2, mesh_2)

    #test A matrices
    @test isapprox(system_1_matrices.A, [-2.39124058 -0.19098 -0.190981; -0.19097978 -0.5 0.94562; 0.190980452 0.945621 -0.5], atol = 0.001)
    @test isapprox(system_2_matrices.A, [-0.68060082 0.0 0.0; 0.0 -0.5 -0.5903; 0.0 -0.5903 -0.5], atol = 0.001)

    #test right hand side vectors
    @test isapprox(system_1_matrices.b, [0.894423 1.18672e-6; 0.447215 -0.89442658; -0.44721 -0.89442776], atol = 0.001)
    @test isapprox(system_2_matrices.b, [0.894423 1.18672e-6; 0.447215 -0.89442658; -0.44721 -0.89442776], atol = 0.001)

    #test A matrices
    @test isapprox(system_1_matrices.A, [-1.1087323 0.19098 -0.190981; 0.19097978 -0.5 0.054366159; -0.190980452 0.054366159 -0.5], atol = 0.001)
    @test isapprox(system_2_matrices.A, [0.18066282145473667 0.0 0.0; 0.0 -0.5 -0.5903; 0.0 -0.5903 -0.5], atol = 0.001)

end

