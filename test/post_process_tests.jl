@testset "Post Processing Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    x, z = x, z = FLOWFoil.naca4()
    coordinates = [x z]
    polar = solve(coordinates, [0.0; 5.0], PlanarProblem(Vortex(Linear()), Dirichlet()))

    Xfoil.set_coordinates(reverse(x), reverse(z))
    cl0, cm0 = Xfoil.solve_alpha(0.0)
    cl5, cm5 = Xfoil.solve_alpha(5.0)

    @test isapprox(polar.lift[1], cl0, atol=1e-5)
    @test isapprox(polar.moment[1], cm0, atol=1e-3)
    @test isapprox(polar.lift[2], cl5, atol=1e-5)
    @test isapprox(polar.moment[2], cm5, atol=1e-3)

    #---------------------------------#
    #           AXISYMMETRIC          #
    #---------------------------------#
    # - Very Basic Test - #

    x, z = FLOWFoil.naca4()
    coordinates = [x z .+ 1.0]
    method = AxisymmetricProblem(Vortex(Constant()), Dirichlet(),[false])

    # Generate Problem Object
    problem = define_problem(method, coordinates, [0.0], [-1.0], [-1.0])

    # Generate Panel Geometry
    panels = generate_panels(method, coordinates)

    # Generate Influence Mesh
    mesh = generate_mesh(method, panels)

    # Assemble Linear System
    system = generate_inviscid_system(method, panels, mesh)

    # Solve Linear System
    solution = solve(system)

    # Post Process Solution
    polar = post_process(method, problem, panels, mesh, solution)

    #TODO: Need to figure out how to test this...
end
