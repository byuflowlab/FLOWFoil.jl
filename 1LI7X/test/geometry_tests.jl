# Test Linear Meshing
@testset "Mesh" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    mesh = FLOWFoil.generate_mesh(x, y)

    xr = [1.0; 0; 0; 1; 1]
    yr = [0.0; 0; 1; 1; 0]
    meshr = FLOWFoil.generate_mesh(xr, yr)

    xd = x
    yd = [-0.01; -0.25; 0; 0.25; 0.01]
    meshd = FLOWFoil.generate_mesh(xd, yd)

    xdr = [0.0; 0; 2]
    ydr = [-1.0; 2; 1]
    meshdr = FLOWFoil.generate_mesh(xdr, ydr)

    # Nodes:
    @testset "Nodes" begin
        nodes = mesh.nodes

        @test isapprox(nodes[1], [1 0])
        @test isapprox(nodes[2], [0.5 -0.5])
        @test isapprox(nodes[3], [0 0])
        @test isapprox(nodes[4], [0.5 0.5])
        @test isapprox(nodes[5], [1 0])
    end
end

# Test distance functions
@testset "Distance Calculations" begin
    x = [1; 0.5; 0; 0.5; 1]
    y = [0; -0.5; 0; 0.5; 0]
    mesh = FLOWFoil.generate_mesh(x, y)
    nodes = mesh.nodes

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
    mesh = FLOWFoil.generate_mesh(x, y)
    nodes = mesh.nodes

    @testset "Vortex Coefficients on Panel" begin
        dmag = sqrt(2) / 2
        r2, r2mag = FLOWFoil.get_r(nodes[2], nodes[1])

        c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[1])

        @test isapprox(
            c1,
            (dmag * log(r2mag) - dmag) / (2 * pi) -
            1.0 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2),
        )
        @test isapprox(c2, 1 / (4 * pi * dmag) * (r2mag^2 * log(r2mag) - 0.5 * r2mag^2))

        r1, r1mag = FLOWFoil.get_r(nodes[1], nodes[2])
        a = sqrt(2) / 2
        c1, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[2])

        @test isapprox(
            c1,
            (1 - a / dmag) * ((a * log(r1mag) - dmag) / (2 * pi)) -
            1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5 * r1mag^2),
        )
        @test isapprox(
            c2,
            (a) * ((a * log(r1mag) - dmag) / (2 * pi * dmag)) +
            1 / (4 * pi * dmag) * (-r1mag^2 * log(r1mag) + 0.5 * r1mag^2),
        )
    end

    @testset "Vortex Coefficients Matrix" begin
        x = [1; 0.5; 0; 0.5; 1]
        y = [0; -0.5; 0; 0.5; 0]
        mesh = FLOWFoil.generate_mesh(x, y)
        nodes = mesh.nodes

        # test psitilde13 as second part of a13
        s22 = sqrt(2) / 2
        h = -s22
        dmag = s22
        ln1 = log(s22)
        ln2 = log(1.0)
        r1mag = s22
        r2mag = 1.0
        theta1 = -pi / 2
        theta2 = -3.0 * pi / 4

        psibar13 = 1.0 / (2 * pi) * (h * (theta2 - theta1) - dmag + dmag * ln2)
        psitilde13 =
            1 / (4 * pi * dmag) *
            (r2mag^2 * ln2 - r1mag^2 * ln1 - 0.5 * r2mag^2 + 0.5 * r1mag^2)

        c1, _ = FLOWFoil.get_vortex_influence(nodes[2], nodes[3], nodes[1])
        @test c1 == psibar13 - psitilde13

        _, c2 = FLOWFoil.get_vortex_influence(nodes[1], nodes[2], nodes[1])
        # Test that coefficient matrix is assembled correctly
        a = FLOWFoil.assemble_vortex_coefficients(mesh,mesh,false)
        @test a[1, 2] == c2 + c1
    end
end
