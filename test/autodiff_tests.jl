@testset "Automatic Differentiation Tests" begin

    # Run some derivatives through each solver

    ######################################################################
    #                                                                    #
    #                               PLANAR                               #
    #                                                                    #
    ######################################################################

    function planar_temp(x0)
        x, z = FLOWFoil.naca4(x0[1], x0[2], x0[3])
        coordinates = [x z]
        polar = solve(coordinates, [0.0; 5.0], PlanarProblem(Vortex(Linear()), Dirichlet()))
        return polar.lift
    end

    x0 = [2.0, 4.0, 12.0]

    # - ForwardDiff - #
    J = ForwardDiff.jacobian(planar_temp, x0)
    @test !any(isnan.(J))

    # - ReverseDiff - #
    # const f_tape = GradientTape(planar_temp, x0)

    # # compile `f_tape` into a more optimized representation
    # const compiled_f_tape = compile(f_tape)

    ######################################################################
    #                                                                    #
    #                            AXISYMMETRIC                            #
    #                                                                    #
    ######################################################################

    function axisym_temp(x0)
        x, z = FLOWFoil.naca4(x0[1], x0[2], x0[3])
        coordinates = [x z .+ 1.0]
        method = AxisymmetricProblem(Vortex(Constant()), Neumann(), [false])

        #TODO: consider putting this in the convenience functions file...
        function solve_axisymmetric(method, coordinates)
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

            return polar
        end

        polar = solve_axisymmetric(method, coordinates)

        return polar.surface_velocity[1]
    end

    # - ForwardDiff - #
    J = ForwardDiff.gradient(axisym_temp, x0)
    @test !any(isnan.(J))

    # - ReverseDiff - #

    ######################################################################
    #                                                                    #
    #                              PERIODIC                              #
    #                                                                    #
    ######################################################################

    function periodic_temp(x0)
        x, z = FLOWFoil.naca4(x0[1], x0[2], x0[3])
        coordinates = [x z .+ 1.0]
        method = PeriodicProblem(Vortex(Constant()), Neumann(), [1.0])

        #TODO: consider putting this in the convenience functions file...
        function solve_periodic(method, coordinates)
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

            return polar
        end

        polar = solve_axisymmetric(method, coordinates)

        return polar.surface_velocity[1]
    end

    # - ForwardDiff - #
    J = ForwardDiff.gradient(axisym_temp, x0)
    @test !any(isnan.(J))

    # - ReverseDiff - #

end
