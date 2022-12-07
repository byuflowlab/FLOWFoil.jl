#=

Convenience functions wrapping problem, system, solution, and post processing steps into single functions for user convenience

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                           SOLVE WRAPPERS                           #
#                                                                    #
######################################################################

#---------------------------------#
#             DEFAULTS            #
#---------------------------------#
#=
    Most users will probably be wanting to do Xfoil-like 2D analyses,
    so we'll have the defaults be xfoil-like implementaions.
=#

"""
    solve(coordinates, angle_of_attack, reynolds, mach, method)

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

**Arguments:**
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `angle_of_attack::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)
- `method::ProblemType` : Type of problem to solve (Planar, axisymmetric, or periodic). Defaults to XFoil-like method: PlanarProblem(Vortex(Linear(),Dirichlet())

**Returns:**
- `polar::Polar` : Polar object according to ProblemType

---

## Other Implementations:
`solve(coordinates, angle_of_attack, reynolds; method)` : assumes mach = 0.0
`solve(coordinates, angle_of_attack; method)` : assumes inviscid and mach = 0.0
`solve(coordinates; method)` : assumes inviscid and mach and angle of attack are 0.0

"""
function solve(
    coordinates,
    angle_of_attack,
    reynolds,
    mach,
    method=PlanarProblem(Vortex(Linear()), Dirichlet()),
)

    # Generate Problem Object
    problem = define_problem(method, coordinates, angle_of_attack, reynolds, mach)

    # Generate Panel Geometry
    panels = generate_panels(method, coordinates)

    # Generate Influence Mesh
    mesh, TEmesh = generate_mesh(method, panels)

    # Assemble Linear System
    system = generate_inviscid_system(method, panels, mesh, TEmesh)

    # Solve Linear System
    solution = solve(system)

    # Post Process Solution
    polar = post_process(method, problem, panels, mesh, solution)

    # Return
    return polar
end

# - Mach = 0.0 - #
function solve(
    coordinates,
    angle_of_attack,
    reynolds,
    method=PlanarProblem(Vortex(Linear()), Dirichlet()),
) end

# - Inviscid; Mach = 0.0 - #
function solve(
    coordinates, angle_of_attack, method=PlanarProblem(Vortex(Linear()), Dirichlet())
)
    return solve(coordinates, angle_of_attack, [-1.0], [-1.0], method)
end

# - Inviscid; AoA = Mach = 0.0 - #
function solve(coordinates, method=PlanarProblem(Vortex(Linear()), Dirichlet()))
    return solve(coordinates, [0.0], [-1.0], [-1.0], method)
end

# TODO: consider also adding convenience functions for the other implementations.
