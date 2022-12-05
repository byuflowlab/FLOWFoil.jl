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
"""
    solve(coordinates, angle_of_attack, reynolds, mach, problemtype)

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

**Arguments:**
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `angle_of_attack::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)
- `problemtype::ProblemType` : Type of problem to solve (planar, axisymmetric, or periodic).

**Returns:**
- `polar::Polar` : Polar object according to ProblemType

---

## Other Implementations:
`solve(coordinates, angle_of_attack, reynolds, problemtype)` : assumes mach = 0.0
`solve(coordinates, angle_of_attack, problemtype)` : assumes inviscid and mach = 0.0
`solve(coordinates, problemtype)` : assumes inviscid and mach and angle of attack are 0.0

"""
function solve(coordinates, angle_of_attack, reynolds, mach, pt::ProblemType)

    # Generate Problem Object
    problem = define_problem(pt, coordinates, angle_of_attack, reynolds, mach)

    # Generate Panel Geometry
    panels = generate_panels(pt, coordinates)

    # Generate Influence Mesh
    mesh, TEmesh = generate_mesh(pt, panels)

    # Assemble Linear System
    system = generate_inviscid_system(pt, panels, mesh, TEmesh)

    # Solve Linear System
    return solution = solve(system)

    # # Post Process Solution
    # polar = post_process(problem, panels, mesh, solution)

    # # Return
    # return polar
    return nothing
end

# - Mach = 0.0 - #
function solve(coordinates, angle_of_attack, reynolds, pt::ProblemType)
    return solve(coordinates, angle_of_attack, reynolds, [-1.0], pt)
end

# - Inviscid; Mach = 0.0 - #
function solve(coordinates, angle_of_attack, pt::ProblemType)
    return solve(coordinates, angle_of_attack, [-1.0], [-1.0], pt)
end

# - Inviscid; AoA = Mach = 0.0 - #
function solve(coordinates, pt::ProblemType)
    return solve(coordinates, [0.0], [-1.0], [-1.0], pt)
end
