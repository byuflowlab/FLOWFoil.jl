#=
Solve Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                          INVISCID SOLVER                           #
#                                                                    #
######################################################################

"""
    solve(problem)

Solve problem defined by the input Problem object and return the solution in a Solution object.

**Arguments:**
- `problem::Problem` : Problem to solve

**Returns:**
 - `solution::{InviscidSolution or ViscousSolution}` : returns solution of type matching viscous flag in problem.
"""
function solve(problem)

    # Check viscosity and solve accordingly
    if problem.reynolds != nothing && problem.reynolds > 0.0
        @error "No viscous implementation."
    else
        solution = solve_inviscid(problem)
    end

    return solution
end

"""
    InviscidSolution{TM,TF,TD}

**Fields:**
 - `x::Array{Float}` : solution value(s) at each airfoil node.
 - `mesh::Mesh` : Mesh object describing airfoil nodes etc.
 - `system::InviscidSystem` : system object.
"""
struct InviscidSolution{TM,TF,TI,TD}
    x::TF
    mesh::TM
    system::TD
end

"""
    solve_inviscid(problem)

Solves the inviscid problem.

**Arguments:**
 - `problem::Problem` : Problem to solve.  viscous field must be set to false.

**Returns:**
 - `solution::InviscidSolution`
"""
function solve_inviscid(problem)

    # Check to make sure you want the invsicid solution:
    if problem.reynolds != nothing && problem.reynolds > 0.0
        @warn "Viscous mismatch, please set problem.viscous=false if you would like to solve the inviscid system alone."
    end

    # get inviscid system
    inviscid_system = generate_inviscid_system(problem.mesh, problem.type)

    solution = solve_inviscid_system(inviscid_system, problem.mesh)

    return solution
end

"""
    solve_inviscid(inviscid_system, mesh)

Solve the inviscid_system for the vortex and streamfunction strengths.

Outputs the InviscidSolution object which contains the inviscid_system.

**Arguments:**
- `inviscid_system::InviscidSystem` : inviscid_system to solve.
- `mesh::Mesh` : Mesh defining geometry (to put into solution object)

**Returns:**
 - `solution::InviscidSolution`

"""
function solve_inviscid(inviscid_system, meshes)

    # Solve System
    x = inviscid_system.vcoeffmat \ inviscid_system.bccoeffvec

    return InviscidSolution(x, meshes, inviscid_system)
end

# # MOVE TO POST PROCESSING
#     # Separate Outputs
#     if axisymmetric
#         nk = countkutta(meshes)
#         panelgammas = gammas[1:(end - nk), :]
#         bodystrength = gammas[(end - nk):end, :]
#     else
#         panelgammas = gammas[1:(end - length(inviscid_system.Ns)), :]
#         bodystrength = gammas[(end - (length(inviscid_system.Ns) - 1)):end, :]
#     end
