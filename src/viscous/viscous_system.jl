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
    rSu = parameters.rSu
    gamma_air = parameters.gamma_air
    machinf = problem.machinf
    reynolds = problem.reynolds
    chord = inviscid_solution.mesh.chord

    #if machinf is non-zero, calculate with Karnam-Tsien compressibility correction
    if machinf > 0.0
        # Karman-Tsien beta
        KTbeta = sqrt(1.0 - machinf^2)

        # Karman-Tsien lambda
        KTlambda = machinf^2 / (1.0 + KTbeta)^2

        # Constant Stagnation Enthalpy
        H0 =
            vinf^2 / ((gamma_air - 1.0) * machinf^2) *
            (1.0 + 0.5 * (gamma_air - 1.0) * machinf^2)

        # Temperature Ratio
        TT0 = 1.0 - 0.5 * vinf^2 / H0

        # Sutherland's Law
        SuTT0 = TT0^1.5 * (1.0 + rSu) / (TT0 + rSu)

        sonic_cp =
            2.0 / (gamma_air * machinf^2) * (
                (
                    (1.0 + 0.5 * (gamma_air - 1.0) * machinf^2) /
                    (1.0 + 0.5 * (gamma_air - 1.0))
                )^(gamma_air / (gamma_air - 1.0)) - 1.0
            ) # CANT FIND THIS IN PAPER... seems to only be used for plotting?
    else
        #else, none of these matters.
        KTbeta = 1.0 #sqrt(1) = 1
        KTlambda = 0.0 #0^2 = 0
        H0 = 0.0
        sonic_cp = 0.0

        # and T/T0 = 1.0
        SuTT0 = 1.0
    end

    # calculate stagnation properties
    # stagnation dynamic viscosity
    mu0 = muinf / SuTT0

    # stagnation density
    rho0 = rhoinf * (1.0 + 0.5 * (gamma_air - 1.0) * machinf^2)^(1.0 / (gamma_air - 1.0))

    # return initialized properties
    return Properties(machinf, KTbeta, KTlambda, H0, sonic_cp, rho0, mu0)
end
