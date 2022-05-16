# Test Linear Meshing
@testset "Mesh" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    mesh = FLOWFoil.generatemesh(x, y)

    xr = [1.0; 0; 0; 1; 1]
    yr = [0.0; 0; 1; 1; 0]
    meshr = FLOWFoil.generatemesh(xr, yr)

    xd = x
    yd = [-0.01; -0.25; 0; 0.25; 0.01]
    meshd = FLOWFoil.generatemesh(xd, yd)

    xdr = [0.0; 0; 2]
    ydr = [-1.0; 2; 1]
    meshdr = FLOWFoil.generatemesh(xdr, ydr)

    # Nodes:
    @testset "Nodes" begin
        nodes = mesh.airfoil_nodes

        @test isapprox(nodes[1], [1 0])
        @test isapprox(nodes[2], [0.5 -0.5])
        @test isapprox(nodes[3], [0 0])
        @test isapprox(nodes[4], [0.5 0.5])
        @test isapprox(nodes[5], [1 0])
    end

    #    # Wake:
    #    @testset "Wake" begin
    #        wnodes = mesh.wake_nodes

    #        @test length(wnodes) == 11
    #        @test isapprox(wnodes[1], [1 0])
    #        @test isapprox(wnodes[end], [2 0])
    #        @test isapprox(wnodes[2], [1.1 0])

    #        wnodesr = meshr.wake_nodes
    #        @test isapprox(wnodesr[1], [1 0])
    #        @test isapprox(wnodesr[end], [1 0] .+ sqrt(2) / 2 * [1 -1])
    #        @test isapprox(wnodesr[2], [1 0] .+ sqrt(2) / 2 * [0.1 -0.1])

    #        wnodesd = meshd.wake_nodes
    #        @test isapprox(wnodes, wnodesd)

    #        wnodesdr = meshdr.wake_nodes
    #        @test isapprox(wnodesdr, wnodesr)
    #    end
end

# Test distance functions
@testset "Distance Calculations" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    mesh = FLOWFoil.generatemesh(x, y)
    nodes = mesh.airfoil_nodes

    @testset "Node on Panel" begin
        # Test r calculation
        r2, magr2 = FLOWFoil.get_r(nodes[2], nodes[1])

        @test r2 == [0.5 0.5]
        @test magr2 == sqrt(2) / 2

        # Test panel length calculation
        d, magd = FLOWFoil.get_d(nodes[1], nodes[2])

        @test magd == sqrt(2) / 2

        r1, magr1 = FLOWFoil.get_r(nodes[1], nodes[1])

        # Test height
        h = FLOWFoil.get_h(r1, d, magd)
        @test h == 0.0

        #test tangent distance
        a = FLOWFoil.get_a(r1, d, magd)
        @test a == 0.0

        # Test angle calculations
        theta1 = FLOWFoil.get_theta(h, a)
        theta2 = FLOWFoil.get_theta(h, a)

        @test isapprox(theta2, pi, atol=1e-7)
    end

    @testset "Node off Panel" begin
        # Test r calculation
        r2, magr2 = FLOWFoil.get_r(nodes[3], nodes[1])

        @test r2 == [1.0 0.0]
        @test magr2 == 1.0

        # Test panel length calculation
        d, magd = FLOWFoil.get_d(nodes[2], nodes[3])

        @test magd == sqrt(2) / 2

        r1, magr1 = FLOWFoil.get_r(nodes[2], nodes[1])
        # Test height
        h = FLOWFoil.get_h(r1, d, magd)

        @test isapprox(h, -sqrt(2) / 2)

        # Test length a
        a = FLOWFoil.get_a(r1, d, magd)

        @test a == 0.0
        # Test angle calculations
        theta1 = FLOWFoil.get_theta(h, a)
        theta2 = FLOWFoil.get_theta(h, a, magd)

        @test isapprox(theta1, -pi / 2)
        @test isapprox(theta2, -3 * pi / 4)
    end
end

# Tests for influence coefficients
@testset "Influence Coefficients" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    mesh = FLOWFoil.generatemesh(x, y)
    nodes = mesh.airfoil_nodes

    @testset "Vortex Coefficients on Panel" begin
        dmag = sqrt(2) / 2
        a = sqrt(2) / 2
        r2, r2mag = FLOWFoil.get_r(nodes[2], nodes[1])

        c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[1])

        @test c1 ==
            -dmag / (2 * pi) - (
            -dmag / (2 * pi) + 1 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2)
        )
        @test c2 ==
            -dmag / (2 * pi) +
              1 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2)

        r1, r1mag = FLOWFoil.get_r(nodes[1], nodes[2])

        c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[2])

        @test c1 ==
            1 / (2 * pi) * (-r1mag + r1mag * log(r1mag)) - (
            1 / (2 * pi) * (-r1mag + r1mag * log(r1mag)) +
            1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5 * r1mag^2)
        )

        @test c2 == (
            1 / (2 * pi) * (-r1mag + r1mag * log(r1mag)) +
            1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5r1mag^2)
        )
    end

    @testset "Vortex Coefficients Matrix" begin
        x = [1; 0.5; 0; 0.5; 1]
        y = [0; -0.5; 0; 0.5; 0]
        mesh = FLOWFoil.generatemesh(x, y)
        nodes = mesh.airfoil_nodes
        meshsystem = FLOWFoil.MeshSystem([mesh], [1.0], [0.0], [[0.0; 0.0]])

        # test psitilde12 as first part of a13
        s22 = sqrt(2) / 2
        psibar12 = -s22 / (2 * pi)
        psitilde12 = psibar12 + 1 / (4 * pi * s22) * (s22^2 * log(s22) - 0.5 * s22^2)
        c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[1])
        @test c2 == psitilde12

        # test psibar13 and psitilde13 as second part of a13
        psibar13 = 1 / (2 * pi) * (s22 * (-pi / 4) - s22 + s22 * log(s22))
        psitilde13 =
            psibar13 +
            1 / (4 * pi * s22) * (1.0 * log(1.0) - s22^2 * log(s22) - 0.5 + 0.5 * s22^2)
        c1, c2 = FLOWFoil.get_vortex_influence(nodes[2], nodes[3], nodes[1])
        @test c1 == psibar13 - psitilde13

        # Test that coefficient matrix is assembled correctly
        a = FLOWFoil.assemblematrixa(meshsystem)
        @test a[1, 2] == psitilde12 + psibar13 - psitilde13
    end
end
