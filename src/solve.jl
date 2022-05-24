#=
Solve Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
"""
function solve_inviscid_system(amat, psi_inf, alpha)
    gammas = amat \ psi_inf
    gamma0 = gammas[1:(end - 1), 1]
    gamma90 = gammas[1:(end - 1), 2]
    gammatot = gammas[:, 1] .* cosd(alpha) .+ gammas[:, 2] .* sind(alpha)

    return gamma0, gamma90, gammatot
end

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

"""
"""
function solve_inviscid(problem)

    # Unpack problem components for convenience.
    # mesh system
    ms = problem.meshsystem

    # freestream
    fs = problem.freestream

    # get inviscid system
    gammas, A, rhs = get_inviscid_system(ms)

    # initalize solution
    #TODO: decide format of solution objects

    for i in 1:length(fs.alpha)
        # get solution and post processing for each angle of attack.
        #TODO: obviously...
    end

    return solution
end

"""
"""
function initalize_viscous(problem)

    # solve inviscid problem first TODO
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

"""
"""
function solve_viscous(problem)

    # initialize viscous solution
    solution = initialize_viscous(problem)

    # solve coupled system
    coupled_solve!(solution)

    # post process viscous solution
    post_process_viscous!(solution)

    return solution
end
