#=
Solve Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
10/22 - axisymmetric solver options available.
=#

"""
    solve(problem)

Solve problem defined by the input Problem object and return the solution in a Solution object.

**Arguments:**
- `problem::Problem` : Problem to solve

**Returns:**
 - `solution::{InviscidSolution or ViscousSolution}` : returns solution of type matching viscous flag in problem.
"""
function solve(problem; parameters=nothing)

    # Check viscosity and solve accordingly
    if problem.viscous
        solution = solve_viscous(problem, parameters)
    else
        solution = solve_inviscid(problem)
    end

    return solution
end

###########################################
#######      INVISCID PROBLEM       #######
###########################################

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
    if problem.viscous
        @warn "Viscous mismatch, please set problem.viscous=false if you would like to solve the inviscid system alone."
    end

    # get inviscid system
    inviscidsystem = get_inviscid_system(problem.meshes; axisymmetric=problem.axisymmetric)

    solution = solve_inviscid_system(inviscidsystem, problem.meshes)

    return solution
end

"""
    solve_inviscid_system(inviscidsystem, mesh)

Solve the InviscidSystem for the vortex and streamfunction strengths.

Outputs the InviscidSolution object which contains the inviscidsystem.

**Arguments:**
- `inviscidsystem::InviscidSystem` : InviscidSystem to solve.
- `mesh::PlanarMesh` : PlanarMesh defining geometry (to put into solution object)

**Returns:**
 - `solution::InviscidSolution`

"""
function solve_inviscid_system(inviscidsystem, meshes)

    # Solve System
    gammas = inviscidsystem.vcoeffmat \ inviscidsystem.bccoeffvec

    # Separate Outputs
    panelgammas = gammas[1:(end - length(inviscidsystem.Ns)), :]
    bodystrength = gammas[(end - (length(inviscidsystem.Ns) - 1)):end, :]

    # Generate Solution Object

    solution = InviscidSolution(
        meshes, panelgammas, bodystrength, inviscidsystem.Ns, inviscidsystem
    )

    return solution
end

############################################
########      VISCOUS  PROBLEM       #######
############################################

#"""
#    solve_viscous(problem)

#Solves the viscous problem.

#**Arguments:**
#- `problem::Problem` : Problem Definition.

#**Returns:**
# - `solution::ViscousSolution`
#"""
#function solve_viscous(problem; parameters=nothing)

#    # Check to make sure you want the invsicid solution:
#    if !problem.viscous
#        @warn(
#            "Viscous mismatch, please set problem.viscous=true if you would like to solve the viscous problem."
#        )
#    end

#    # initialize viscous solution
#    solution = initialize_viscous(problem)

#    # solve coupled system
#    coupled_solve!(solution)

#    return solution
#end

#"""
#    initalize_viscous(problem)

#Initialized viscous solution (solves invscid problem, initialized wake and boundary layer, etc.)

#**Arguments:**
# - problem::Problem` : Problem to solve.

#**Returns:**
# - solution::ViscousSolution` : Initialized ViscousSolution
#"""
#function initalize_viscous(problem; parameters=nothing)

#    # solve inviscid problem first
#    inviscid_solution = solve_inviscid(problem)

#    # initialize parameters if needed
#    if parameters == nothing
#        parameters = defaultparameters(;
#            muinf=inviscid_solution.mesh.chord / problem.reynolds
#        )
#    end

#    # set thermodynamic properties TODO
#    properties = initialize_properties(problem, parameters, inviscid_solution)

#    # locate stagnation point TODO

#    # split into pressure/suction sides TODO

#    # generate wake TODO
#    wake = generate_wake(problem, inviscid_solution)

#    # set dead air space behind blunt trailing edges. TODO

#    # initialize edge velocity TODO

#    # initialize boundary layer TODO

#    # update stagnation point location based on edge velocities of initialized boundary layer. TODO

#    # initialize viscous solution object TODO
#    solution = ViscousSolution()

#    return solution
#end
