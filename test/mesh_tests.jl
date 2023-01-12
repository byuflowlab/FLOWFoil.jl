@testset "Planar Mesh Tests" begin

    # - Very Basic Test - #

    x = [1.0; 0.0; 1.0]
    z = [0.0; 0.0; eps()]
    coordinates = [x z]
    panels = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    mesh, TEmesh = generate_mesh(PlanarProblem(Vortex(Linear()), Dirichlet()), panels)

    @test mesh.panel_indices == [1:2]
    @test mesh.node_indices == [1:3]
    @test mesh.mesh2panel == [1; 2]
    @test mesh.chord == 1.0
    @test mesh.panel_length == [1.0; 1.0]
    @test isapprox(mesh.r1, [0.0 1.0; 1.0 0.0; eps() 1.0])
    @test mesh.lnr1 == zeros(3, 2)
    @test isapprox(mesh.r1normal, [0.0 -eps(); 0.0 0.0; 0.0 0.0])
    @test mesh.r1tangent == [0.0 1.0; 1.0 0.0; 0.0 1.0]
    @test mesh.theta1 == [pi 0.0; 0.0 pi; pi 0.0]
    @test isapprox(mesh.r2, [1.0 eps(); 0.0 1.0; 1.0 0.0])
    @test mesh.lnr2 == zeros(3, 2)
    @test mesh.theta2 == [pi 0.0; 0.0 pi; pi 0.0]

    @test TEmesh.blunt_te == Bool[0]
    @test TEmesh.trailing_edge_gap == [0.0]
    @test TEmesh.tdp == [1.0]
    @test TEmesh.txp == [0.0]

    #---------------------------------#
    #         Multi Body Test         #
    #---------------------------------#
    x1 = zeros(3)
    z1 = [-0.5; 0.0; 0.5]

    x2 = [1.0; 1.5; 2.0]
    z2 = zeros(3)

    coordinates = ([x1 z1], [x2 z2])
    panels = generate_panels(PlanarProblem(Vortex(Linear()), Dirichlet()), coordinates)
    mesh, TEmesh = generate_mesh(PlanarProblem(Vortex(Linear()), Dirichlet()), panels)

    @test mesh.panel_indices == [[1:2]; [3:4]]
    @test mesh.node_indices == [[1:3]; [4:6]]
    @test mesh.mesh2panel == [1; 2; 1; 2]
    @test mesh.chord == 2.0
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

    @test TEmesh.blunt_te == Bool[1; 1]
    @test TEmesh.trailing_edge_gap == [0.0; 0.0]
    @test TEmesh.tdp == [1.0; 1.0]
    @test TEmesh.txp == [0.0; 0.0]
end

@testset "Axisymmetric Mesh Tests" begin

    # - Very Basic Test - #

    x = [0.0; 1.0; 2.0]
    z = [0.0; 1.0; 2.0]
    coordinates = [x z]
    panels = generate_panels(
        AxisymmetricProblem(Vortex(Constant()), Neumann(), [true]), coordinates
    )
    mesh = generate_mesh(AxisymmetricProblem(Vortex(Constant()), Neumann(), [true]), panels)

    @test mesh.nbodies == 1
    @test mesh.panel_indices == [1:2]
    @test mesh.x == [0.0 (0.5 - 1.5)/1.5; 1.0/0.5 0.0]
    @test mesh.r == [1.0 0.5/1.5; 1.5/0.5 1.0]
    @test mesh.m == [1.0 0.6; 0.6 1.0]

    # - Multi Body Test - #

    x1 = [0.0; 0.5; 1.0]
    r1 = [0.0; 0.5; 0.0]

    x2 = [1.0; 0.5; 0.0; 0.5; 1.0]
    r2 = [1.5; 1.0; 1.5; 2.0; 1.5]

    coordinates = ([x1 r1], [x2 r2])
    panel_array = generate_panels(
        AxisymmetricProblem(Vortex(Constant()), Neumann(), [true, false]), coordinates
    )
    mesh = generate_mesh(
        AxisymmetricProblem(Vortex(Constant()), Neumann(), [true, false]), panel_array
    )

    @test mesh.nbodies == 2
    @test mesh.panel_indices == [[1:2]; [3:6]]
    @test mesh.mesh2panel == [1;2;1;2;3;4]
    @test mesh.x == [
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
    ]
    ri = [0.25; 0.25; 1.25; 1.25; 1.75; 1.75]
    @test isapprox(mesh.r, ri * (1.0 ./ ri'))
    x = [
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (0.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(1.0 / 4.0) (-1.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0) (-1.0 / 2.0)/(7.0 / 4.0)
        (1.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(1.0 / 4.0) (0.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(5.0 / 4.0) (1.0 / 2.0)/(7.0 / 4.0) (0.0 / 2.0)/(7.0 / 4.0)
    ]
    r = ri * (1.0 ./ ri')
    @test isapprox(mesh.m , 4.0*r./(x.^2 .+ (r .+ 1.0).^2))
end

@testset "Periodic Mesh Tests" begin

    # - Very Basic Test - #

    x = [0.0; 1.0; 2.0]
    z = [0.0; 1.0; 2.0]
    coordinates = [x z]
    panels = generate_panels(
        PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0), coordinates
    )
    mesh = generate_mesh(PeriodicProblem(Vortex(Constant()), Neumann(), 1.0, 0.0), panels)

    @test mesh.nbodies == 1
    @test mesh.panel_indices == [1:2]
    @test mesh.x == [0.0 -1.0; 1.0 0.0]
    @test mesh.y == [0.0 -1.0; 1.0 0.0]
end
