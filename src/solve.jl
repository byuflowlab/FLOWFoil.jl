#=
Solve Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
"""
function solve(problem)

    # Check viscosity and solve accordingly
    if problem.viscous
        solution = solve_viscous(problem)
    else
        solution = solve_inviscid(problem)
    end

    return solution
end

###########################################
#######      INVISCID PROBLEM       #######
###########################################

"""
"""
function solve_inviscid(problem)

    # Generate Mesh
    mesh = generate_mesh(problem.coordinates)

    # get inviscid system
    inviscidsystem = get_inviscid_system(mesh)

    # solve inviscid system
    solution = solve_inviscid_system(inviscidsystem, mesh; debug=problem.debug)

    return solution
end

"""
"""
function solve_inviscid_system(inviscidsystem, mesh; debug=false)

    # Solve System
    gammas = inviscidsystem.vcoeffmat \ inviscidsystem.bccoeffvec

    # Separate Outputs
    panelgammas = gammas[1:(end - 1), :]
    psi0 = gammas[end, :]

    # Generate Solution Object
    if debug
        solution = InviscidSolution(mesh, panelgammas, psi0, inviscidsystem)
    else
        solution = InviscidSolution(mesh, panelgammas, psi0, nothing)
    end

    return solution
end

###########################################
#######      VISCOUS  PROBLEM       #######
###########################################

"""
"""
function solve_viscous(problem)

    # initialize viscous solution
    solution = initialize_viscous(problem)

    # solve coupled system
    coupled_solve!(solution)

    return solution
end

"""
"""
function initalize_viscous(problem)

    # solve inviscid problem first
    inviscid_solution = solve_inviscid(problem)

    # set thermodynamic properties TODO
    properties = set_properties(problem)

    # generate wake TODO
    wake = generate_wake(problem, inviscid_solution)

    # locate stagnation point TODO

    # split into pressure/suction sides TODO

    # set dead air space behind blunt trailing edges. TODO

    # initialize edge velocity TODO

    # initialize boundary layer TODO

    # update stagnation point location based on edge velocities of initialized boundary layer. TODO

    # initialize viscous solution object TODO
    solution = ViscousSolution()

    return solution
end
