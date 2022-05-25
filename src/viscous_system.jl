#=
Boudary Layer Integral Method with Transition Models

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    initialize_properties(problem, parameters, inviscid_solution)

Initialize thermodynamic properties for viscous solver.

**Arguments:**
 - 'problem::Problem' : Problem to solve.
 - 'parameters::Parameters' : Solver parameters.
 - 'inviscid_solution::InviscidSolution' : Solution to inviscid problem.

**Returns:**
 - 'properties::Properties' : Thermodynamic properties object
"""
function initialize_properties(problem, parameters, inviscid_solution)

    #rename for convenience
    rhoinf = parameters.rhoinf
    vinf = parameters.vinf
    muinf = parameters.muinf
    Tsrat = parameters.Tsrat
    machinf = problem.machinf
    reynolds = problem.reynolds
    chord = inviscid_solution.mesh.chord

    # Get gamma magnitude at nodes
    mag_gammas = get_vortex_magnitudes(inviscid_solution, problem.angleofattack)
    mag_gammas_m1 = mag_gammas .- 1.0

    #if machinf is non-zero, apply compressibility items
    if machinf > 0.0
        KTb = sqrt(1.0 - machinf^2)
        KTl = machinf^2 / (1.0 + KTb)^2
        H0 =
            (1.0 .+ 0.5 * mag_gammas_m1 * machinf^2) * vinf^2.0 ./
            (mag_gammas_m1 * machinf^2) #TODO: something is weird here. need to look up matlab order of operations and see which operators are elementwise by default. This should result in a single value

        Tr = 1.0 - 0.5 * vinf^2 / H0
        finf = Tr^1.5 * (1.0 + Tsrat) / (Tr + Tsrat)
        #cps = 2/(g*Minf^2)*(((1+0.5*gmi*Minf^2)/(1+0.5*gmi))^(g/gmi) - 1) #is this needed?
    else
        KTb = 0.0
        KTl = 0.0
        H0 = 0.0
        finf = 1.0
    end

    #calculate stagnation properties
    mu0 = muinf / finf
    rho0 = rhoinf * (1.0 .+ 0.5 * mag_gammas_m1 * machinf^2) .^ (1.0 ./ mag_gammas_m1) #TODO: similar: need to look up matlab stuff, this should be a single value in the end.  Which operators are elementwise and which arent?

    # return initialized properties
    return Properties(machinf, KTb, KTl, H0, rho0, mu0)
end
