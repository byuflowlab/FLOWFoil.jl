@testset "Planar Panel Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#

    # - Very Basic Test - #
    x = [1.0; 2.0; 3.0; 4.0]
    z = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x z]
    panel_array = FLOWFoil.generate_panel_geometry(Mfoil(), coordinates)
    panels = panel_array
    @test panels.npanels == 3
    @test panels.panel_edges[1, :, :] == [1.0 0.0; 2.0 0.0]
    @test panels.panel_edges[2, :, :] == [2.0 0.0; 3.0 0.0]
    @test panels.panel_edges[3, :, :] == [3.0 0.0; 4.0 0.0]
    @test panels.panel_lengths == ones(3)
    @test panels.panel_vectors == [1.0 0.0; 1.0 0.0; 1.0 0.0]

    # - Multi-body Tests - #
    x1 = [1.0; 0.5; 0.0; 0.5; 1.0]
    z1 = [0.0; -0.5; 0.0; 0.5; 0.0]

    x2 = 1.25 .+ [1.0; 0.5; 0.0; 0.5; 1.0]
    z2 = [0.0; -0.5; 0.0; 0.5; 0.0]

    coordinates = ([x1 z1], [x2 z2])
    panel_array = FLOWFoil.generate_panel_geometry(Mfoil(), coordinates)

    @test panel_array[1].npanels == 4
    @test panel_array[1].panel_edges[1, :, :] == [1.0 0.0; 0.5 -0.5]
    @test panel_array[1].panel_edges[2, :, :] == [0.5 -0.5; 0.0 0.0]
    @test panel_array[1].panel_edges[3, :, :] == [0.0 0.0; 0.5 0.5]
    @test panel_array[1].panel_edges[4, :, :] == [0.5 0.5; 1.0 0.0]
    @test panel_array[1].panel_lengths == sqrt(0.5) * ones(4)
    @test panel_array[1].panel_vectors == [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]

    @test panel_array[2].npanels == 4
    @test panel_array[2].panel_edges[1, :, :] == [1.25 0.0] .+ [1.0 0.0; 0.5 -0.5]
    @test panel_array[2].panel_edges[2, :, :] == [1.25 0.0] .+ [0.5 -0.5; 0.0 0.0]
    @test panel_array[2].panel_edges[3, :, :] == [1.25 0.0] .+ [0.0 0.0; 0.5 0.5]
    @test panel_array[2].panel_edges[4, :, :] == [1.25 0.0] .+ [0.5 0.5; 1.0 0.0]
    @test panel_array[2].panel_lengths == sqrt(0.5) * ones(4)
    @test panel_array[2].panel_vectors == [-0.5 -0.5; -0.5 0.5; 0.5 0.5; 0.5 -0.5]
end


