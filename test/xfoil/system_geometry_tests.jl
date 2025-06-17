# Test Linear Meshing
@testset "Mesh" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [x y])
    mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    xr = [1.0; 0; 0; 1; 1]
    yr = [0.0; 0; 1; 1; 0]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [xr yr])
    meshr = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    xd = x
    yd = [-0.01; -0.25; 0; 0.25; 0.01]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [xd yd])
    meshd = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    xdr = [0.0; 0; 2]
    ydr = [-1.0; 2; 1]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [xdr ydr])
    meshdr = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    # Nodes:
    @testset "Nodes" begin
        nodes = mesh.nodes

        @test isapprox(nodes[1, :], [1; 0])
        @test isapprox(nodes[2, :], [0.5; -0.5])
        @test isapprox(nodes[3, :], [0; 0])
        @test isapprox(nodes[4, :], [0.5; 0.5])
        @test isapprox(nodes[5, :], [1; 0])
    end
end

# Test distance functions
@testset "Distance Calculations" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [x y])
    mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)
    nodes = mesh.nodes

    @testset "Node on Panel" begin
        # Test r calculation
        r2, magr2 = FLOWFoil.get_r(nodes[2, :], nodes[1, :])

        @test r2 == [0.5; 0.5]
        @test magr2 == sqrt(2) / 2

        # Test panel length calculation
        d, magd = FLOWFoil.get_d(nodes[1, :], nodes[2, :])

        @test magd == sqrt(2) / 2

        r1, magr1 = FLOWFoil.get_r(nodes[1, :], nodes[1, :])

        # TODO: figure out what these got renamed to
        # # Test height
        # h = FLOWFoil.get_h(r1, d, magd)
        # @test h == 0.0
        #
        # #test tangent distance
        # a = FLOWFoil.get_a(r1, d, magd)
        # @test a == 0.0
        #
        # # Test angle calculations
        # theta1 = FLOWFoil.get_theta(h, a)
        # theta2 = FLOWFoil.get_theta(h, a)
        #
        # @test isapprox(theta2, pi, atol=1e-7)
    end

    @testset "Node off Panel" begin
        # Test r calculation
        r2, magr2 = FLOWFoil.get_r(nodes[3, :], nodes[1, :])

        @test r2 == [1.0; 0.0]
        @test magr2 == 1.0

        # Test panel length calculation
        d, magd = FLOWFoil.get_d(nodes[2, :], nodes[3, :])

        @test magd == sqrt(2) / 2

        r1, magr1 = FLOWFoil.get_r(nodes[2, :], nodes[1, :])

        # TODO: figure out what these got renamed to
        # # Test height
        # h = FLOWFoil.get_h(r1, d, magd)
        # @test isapprox(h, -sqrt(2) / 2)
        #
        # # Test length a
        # a = FLOWFoil.get_a(r1, d, magd)
        #
        # @test a == 0.0
        # # Test angle calculations
        # theta1 = FLOWFoil.get_theta(h, a)
        # theta2 = FLOWFoil.get_theta(h, a, magd)
        #
        # @test isapprox(theta1, -pi / 2)
        # @test isapprox(theta2, -3 * pi / 4)
    end
end

#TODO: update these to new functions
# # Tests for influence coefficients
# @testset "Influence Coefficients" begin
#     x = [1; 0.5; 0; 0.5; 1]
#     y = [0; -0.5; 0; 0.5; 0]
#     panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [x y])
#     mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)
#     nodes = mesh.nodes
#
#     @testset "Vortex Coefficients on Panel" begin
#         dmag = sqrt(2) / 2
#         r2, r2mag = FLOWFoil.get_r(nodes[2, :], nodes[1, :])
#
#         c1, c2 = FLOWFoil.get_vortex_influence(nodes[1, :], nodes[2, :], nodes[1, :])
#
#         @test isapprox(
#             c1,
#             (dmag * log(r2mag) - dmag) / (2 * pi) -
#             1.0 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2),
#         )
#         @test isapprox(c2, 1 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2))
#
#         r1, r1mag = FLOWFoil.get_r(nodes[1], nodes[2])
#         a = sqrt(2) / 2
#         c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[2])
#
#         @test isapprox(
#             c1,
#             (1 - a / dmag) * ((a * log(r1mag) - dmag) / (2 * pi)) -
#             1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5 * r1mag^2),
#         )
#         @test isapprox(
#             c2,
#             (a) * ((a * log(r1mag) - dmag) / (2 * pi * dmag)) +
#             1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5 * r1mag^2),
#         )
#     end
#
#     @testset "Vortex Coefficients Matrix" begin
#         x = [1; 0.5; 0; 0.5; 1]
#         y = [0; -0.5; 0; 0.5; 0]
#         panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), [x y])
#         mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)
#         nodes = mesh.nodes
#
#         # test psitilde13 as second part of a13
#         s22 = sqrt(2) / 2
#         h = -s22
#         dmag = s22
#         ln1 = log(s22)
#         ln2 = log(1.0)
#         r1mag = s22
#         r2mag = 1.0
#         theta1 = -pi / 2
#         theta2 = -3.0 * pi / 4
#
#         psibar13 = 1.0 / (2 * pi) * (h * (theta2 - theta1) - dmag + dmag * ln2)
#         psitilde13 =
#             1 / (4 * pi * dmag) *
#             (r2mag^2 * ln2 - r1mag^2 * ln1 - 0.5 * r2mag^2 + 0.5 * r1mag^2)
#
#         c1, _ = FLOWFoil.get_vortex_influence(nodes[2], nodes[3], nodes[1])
#         @test c1 == psibar13 - psitilde13
#
#         _, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[1])
#         # Test that coefficient matrix is assembled correctly
#         a = FLOWFoil.assemble_vortex_coefficients(mesh, mesh, false)
#         @test a[1, 2] == c2 + c1
#     end
# end

@testset "Planar Mesh Tests" begin

    # - Very Basic Test - #

    x = [1.0; 0.0; 1.0]
    y = [0.0; 0.0; eps()]
    coordinates = [x y]
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), coordinates)
    mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    @test mesh.panel_indices == [1:2]
    @test mesh.node_indices == [1:3]
    @test mesh.mesh2panel == [1; 2]
    @test mesh.chord_length == 1.0
    @test mesh.panel_length == [1.0; 1.0]
    @test isapprox(mesh.r1, [0.0 1.0; 1.0 0.0; eps() 1.0])
    @test mesh.lnr1 == zeros(3, 2)
    @test isapprox(mesh.r1normal, [0.0 -eps(); 0.0 0.0; 0.0 0.0])
    @test mesh.r1tangent == [0.0 1.0; 1.0 0.0; 0.0 1.0]
    @test mesh.theta1 == [pi 0.0; 0.0 pi; pi 0.0]
    @test isapprox(mesh.r2, [1.0 eps(); 0.0 1.0; 1.0 0.0])
    @test mesh.lnr2 == zeros(3, 2)
    @test mesh.theta2 == [pi 0.0; 0.0 pi; pi 0.0]

    @test mesh.TE_geometry.blunt_te == Bool[0]
    @test mesh.TE_geometry.panel_length == [0.0]
    @test mesh.TE_geometry.tdp == [1.0]
    @test mesh.TE_geometry.txp == [0.0]

    #---------------------------------#
    #         Multi Body Test         #
    #---------------------------------#
    x1 = zeros(3)
    z1 = [-0.5; 0.0; 0.5]

    x2 = [1.0; 1.5; 2.0]
    z2 = zeros(3)

    coordinates = ([x1 z1], [x2 z2])
    panel_geometry = FLOWFoil.generate_panel_geometry(Xfoil(), coordinates)
    mesh = FLOWFoil.generate_system_geometry(Xfoil(), panel_geometry)

    @test mesh.panel_indices == [[1:2]; [3:4]]
    @test mesh.node_indices == [[1:3]; [4:6]]
    @test mesh.mesh2panel == [1; 2; 1; 2]
    @test mesh.chord_length == 2.0
    @test all(mesh.panel_length .== 0.5 .* ones(4))
    @test all(
        mesh.r1 .== [
            0.0 0.5 sqrt(1.25) sqrt(2.5)
            0.5 0.0 1.0 1.5
            1.0 0.5 sqrt(1.25) sqrt(2.5)
            sqrt(1.25) 1.0 0.0 0.5
            sqrt(2.5) 1.5 0.5 0.0
            sqrt(4.25) 2.0 1.0 0.5
        ],
    )
    @test all(
        mesh.lnr1 .== [
            0.0 log(0.5) log(sqrt(1.25)) log(sqrt(2.5))
            log(0.5) 0.0 log(1.0) log(1.5)
            log(1.0) log(0.5) log(sqrt(1.25)) log(sqrt(2.5))
            log(sqrt(1.25)) log(1.0) 0.0 log(0.5)
            log(sqrt(2.5)) log(1.5) log(0.5) 0.0
            log(sqrt(4.25)) log(2.0) log(1.0) log(0.5)
        ],
    )
    @test all(
        mesh.r1normal .== [
            0.0 0.0 -0.5 -0.5
            0.0 0.0 0.0 0.0
            0.0 0.0 0.5 0.5
            -1.0 -1.0 0.0 0.0
            -1.5 -1.5 0.0 0.0
            -2.0 -2.0 0.0 0.0
        ],
    )
    @test all(
        mesh.r1tangent .== [
            0.0 -0.5 -1.0 -1.5
            0.5 0.0 -1.0 -1.5
            1.0 0.5 -1.0 -1.5
            0.5 0.0 0.0 -0.5
            0.5 0.0 0.5 0.0
            0.5 0.0 1.0 0.5
        ],
    )
    @test all(
        mesh.theta1 .== [
            pi -pi -(pi / 2.0 + atan(1.0, 0.5)) -(pi / 2.0 + atan(1.5, 0.5))
            0.0 pi pi pi
            0.0 0.0 pi / 2.0+atan(1.0, 0.5) pi / 2.0+atan(1.5, 0.5)
            -atan(1.0, 0.5) -pi/2 pi pi
            -atan(1.5, 0.5) -pi/2 0.0 pi
            -atan(2.0, 0.5) -pi/2 0.0 0.0
        ],
    )

    @test all(
        mesh.r2 .== [
            0.5 1.0 sqrt(2.5) sqrt(4.25)
            0.0 0.5 1.5 2.0
            0.5 0.0 sqrt(2.5) sqrt(4.25)
            1.0 sqrt(1.25) 0.5 1.0
            1.5 sqrt(2.5) 0.0 0.5
            2.0 sqrt(4.25) 0.5 0.0
        ],
    )
    @test all(
        mesh.lnr2 .== [
            log(0.5) log(1.0) log(sqrt(2.5)) log(sqrt(4.25))
            0.0 log(0.5) log(1.5) log(2.0)
            log(0.5) 0.0 log(sqrt(2.5)) log(sqrt(4.25))
            log(1.0) log(sqrt(1.25)) log(0.5) log(1.0)
            log(1.5) log(sqrt(2.5)) 0.0 log(0.5)
            log(2.0) log(sqrt(4.25)) log(0.5) 0.0
        ],
    )
    @test all(
        isapprox.(
            mesh.theta2,
            [
                pi -pi -(pi / 2.0 + atan(1.5, 0.5)) -(pi / 2.0 + atan(2.0, 0.5))
                0.0 pi pi pi
                0.0 0.0 pi / 2.0+atan(1.5, 0.5) pi / 2.0+atan(2.0, 0.5)
                -pi/2 -(pi / 2 + atan(0.5, 1.0)) pi pi
                -pi/2 -(pi / 2 + atan(0.5, 1.5)) 0.0 pi
                -pi/2 -(pi / 2 + atan(0.5, 2.0)) 0.0 0.0
            ],
        ),
    )

    @test mesh.TE_geometry.blunt_te == Bool[1; 1]
    @test mesh.TE_geometry.panel_length == [0.0; 0.0]
    @test mesh.TE_geometry.tdp == [1.0; 1.0]
    @test mesh.TE_geometry.txp == [0.0; 0.0]

    #TODO need to put together better trailing edge mesh test with geometry that actually has a proper trailing edge gap.
end

