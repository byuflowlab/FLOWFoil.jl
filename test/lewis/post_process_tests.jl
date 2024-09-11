
@testset "Post Processing Tests" begin

    # - Very Basic Test - #

    x, z = FLOWFoil.naca4()
    coordinates = [x z .+ 1.0]
    method = AxisymmetricProblem(Vortex(Constant()), Dirichlet(), [false])

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
