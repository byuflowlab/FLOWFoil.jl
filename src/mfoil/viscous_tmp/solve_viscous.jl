
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
