#=
Inviscid Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 16 May 2022

Change Log:
=#

"""
"""
function inviscid_post(inviscid_solution, angleofattack)

    # rename fields for convenience.
    gamma0 = inviscid_solution.panelgammas[:, 1]
    gamma90 = inviscid_solution.panelgammas[:, 2]
    nodes = inviscid_solution.mesh.airfoil_nodes
    N = length(gamma0)

    #calculate chord
    chord = maximum(getindex.(nodes, 1)) - minimum(getindex.(nodes, 1))

    # Get Velocity at NODES
    vti = [gamma0[i] * cosd(angleofattack) + gamma90[i] * sind(angleofattack) for i in 1:N]

    # Get Pressure at NODES
    cpi = 1.0 .- vti .^ 2

    # Get Mean Pressure at PANEL MIDPOINTS
    cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0

    # Get Panel Vectors
    dn = nodes[2:end] .- nodes[1:(end - 1)]

    # Get Lift Coefficient
    cl =
        sum([
            cpibar[i] * (-sind(angleofattack) * dn[i][2] - cosd(angleofattack) * dn[i][1])
            for i in 1:(N - 1)
        ]) / chord

    # Get Drag Coefficient
    cdp = 0.0
    cdi =
        sum([
            cpibar[i] * (cosd(angleofattack) * dn[i][2] - sind(angleofattack) * dn[i][1])
            for i in 1:(N - 1)
        ]) / chord
    cd = cdi

    # Get Moment Coefficient
    # leading edge location for moment reference)
    # _, leadingedgeidx = findmin(getindex.(nodes, 1))
    # x0 = nodes[leadingedgeidx][1]
    # z0 = nodes[leadingedgeidx][2]
    #quarter chord location (moment reference location for inviscid case)
    x0 = chord/4.0
    z0 = 0.0

    #define portions of moment calculation
    M = [2 1; 1 2] ./ 6.0
    dxddmi = [
        dn[i][1] * (nodes[i][1] - x0) + dn[i][2] * (nodes[i][2] - z0) for i in 1:(N - 1)
    ]
    dxddmip1 = [
        dn[i][1] * (nodes[i + 1][1] - x0) + dn[i][2] * (nodes[i + 1][2] - z0) for
        i in 1:(N - 1)
    ]

    #calculate moment coefficient about leading edge
    cmi =
        sum([[cpi[i] cpi[i + 1]] * M * [dxddmi[i]; dxddmip1[i]] for i in 1:(N - 1)]) / chord^2

    # Create Polar Object
    polar = Polar(cl, cd, cdp, cdi, cmi[1], vti, cpi)

    return polar
end

"""
    get_vortex_magnitudes(inviscid_solution,angleofattack)

Calculate the vortex strength magnitudes at the airfoil nodes for a given angle of attack.

**Arguments:**
 - 'inviscid_solution::InviscidSolution' : the inviscid solution from which to find the vortex magnitudes.
 - 'angleofattack::Float' : the angle of attack in degrees.
"""
function get_vortex_magnitudes(inviscid_solution, angleofattack)

    #rename for convenience
    gammas = inviscid_solution.panelgammas

    return [
        gammas[i, 1] * cosd(angleofattack) + gammas[i, 2] * sind(angleofattack) for
        i in 1:length(gammas)
    ]
end
