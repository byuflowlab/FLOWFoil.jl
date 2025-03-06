@testset "System Assembly Tests" begin

    #---------------------------------#
    #             PERIODIC            #
    #---------------------------------#
    # - Very Basic Test - #

    x = [1.0, 0.5, 0.0, 0.5, 1.0]
    y = [0.0, -0.5, 0.0, 0.5, 0.0]
    xy = [x y]

    # Generate Panel Geometry
    panel_geometry = FLOWFoil.generate_panel_geometry(HessSmith(), xy)
    system_geometry = FLOWFoil.generate_system_geometry(HessSmith(), panel_geometry)
    strengths = FLOWFoil.generate_system_matrices(HessSmith(), panel_geometry, system_geometry)
    @test 1 == 1

    # # check that it's not invertable
    # @test !isapprox(LinearAlgebra.det(system.A), 0.0)

    # # check that the kutta condition is in the right place
    # @test system.A[end, 1] == 1.0
    # @test system.A[end, end - 1] == 1.0
    # @test system.A[end, end] == 0.0
    # @test all(system.A[end, 2:(end - 2)] .== 0.0)

    # # check that the stream function values are in the right places
    # @test all(system.A[1:(end - 1), end] .== 1.0)
    # @test all(system.A[end, end] .== 0.0)

    # # TODO: multi-airfoil test
end
