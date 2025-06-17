@testset "Axisymmetric Panel Tests" begin
    #---------------------------------#
    #           AXISYMMETRIC          #
    #---------------------------------#
    # Very Basic Test
    x = [1.0; 2.0; 3.0; 4.0]
    y = [0.0; 0.0; 0.0; 0.0]
    coordinates = [x y]
    panels = FLOWFoil.generate_panel_geometry(
        Lewis(; body_of_revolution=[true]), coordinates
    )
    @test panels.npanels == 3
    @test panels.panel_center == [1.5 0.0; 2.5 0.0; 3.5 0.0]
    @test panels.panel_length == ones(3)
    @test panels.panel_normal == [0.0 1.0; 0.0 1.0; 0.0 1.0]
    @test panels.panel_curvature == zeros(3)
    @test panels.panel_angle == zeros(3)

    # Multibody test
    x1 = [0.0; 0.5; 1.0]
    r1 = [0.0; 0.5; 0.0]

    x2 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r2 = [1.5; 1.0; 1.5; 2.0; 1.5]

    coordinates = ([x1 r1], [x2 r2])
    panel_array = FLOWFoil.generate_panel_geometry(
        Lewis(; body_of_revolution=[true, false]), coordinates
    )

    @test panel_array[1].npanels == 2
    @test panel_array[1].panel_center == [0.25 0.25; 0.75 0.25]
    @test panel_array[1].panel_length == [sqrt(1 / 2); sqrt(1 / 2)]
    @test isapprox(panel_array[1].panel_normal, [-sqrt(2)/2 sqrt(2)/2; sqrt(2)/2 sqrt(2)/2])
    @test panel_array[1].panel_curvature == [0.0; 0.0]
    @test panel_array[1].panel_angle == [pi / 4.0; -pi / 4.0]

    @test panel_array[2].npanels == 4
    @test panel_array[2].panel_center == [0.75 1.25; 0.25 1.25; 0.25 1.75; 0.75 1.75]
    @test panel_array[2].panel_length == sqrt(1 / 2) .* ones(4)
    @test isapprox(
        panel_array[2].panel_normal,
        [
            sqrt(2)/2 -sqrt(2)/2
            -sqrt(2)/2 -sqrt(2)/2
            -sqrt(2)/2 sqrt(2)/2
            sqrt(2)/2 sqrt(2)/2
        ],
    )
    @test panel_array[2].panel_curvature == [0.0; -1.0 / 8.0; -1.0 / 8.0; 0.0]
    @test panel_array[2].panel_angle ==
        [5.0 * pi / 4.0; 3.0 * pi / 4.0; pi / 4.0; -pi / 4.0]
end
