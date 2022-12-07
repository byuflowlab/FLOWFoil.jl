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
function solve(system::System)
    solution = solve_inviscid(system)

    return solution
end

"""
    InviscidSolution{TM,TF,TD}

**Fields:**
 - `x::Array{Float}` : solution value(s) at each airfoil node.
 - `mesh::Mesh` : Mesh object describing airfoil nodes etc.
 - `system::InviscidSystem` : system object.
"""
struct InviscidSolution{TF,TS<:System}
    x::TF
    system::TS
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
function solve_inviscid(inviscid_system)

    # Solve System
    x = inviscid_system.A \ inviscid_system.b

    return InviscidSolution(x, inviscid_system)
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
