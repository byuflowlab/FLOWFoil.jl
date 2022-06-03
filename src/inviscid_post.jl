#=
Inviscid Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 16 May 2022

Change Log:
=#

"""
"""
function inviscid_post(inviscid_solution, angleofattack; cascade=false)
    M = length(inviscid_solution.Ns)

    # rename fields for convenience.
    gamma0 = inviscid_solution.panelgammas[:, 1]
    gamma90 = inviscid_solution.panelgammas[:, 2]
    Ns = inviscid_solution.Ns
    N = sum(Ns)
    meshes = inviscid_solution.meshes

    #calculate chord
    if cascade
        chord =
            maximum(getindex.(meshes[1].airfoil_nodes, 1)) -
            minimum(getindex.(meshes[1].airfoil_nodes, 1))
    else
        chord =
            maximum([maximum(getindex.(meshes[i].airfoil_nodes, 1)) for i in 1:M]) - minimum([minimum(getindex.(meshes[i].airfoil_nodes, 1)) for i in 1:M])
    end

    # Get Velocity at NODES
    vti = [gamma0[i] * cosd(angleofattack) + gamma90[i] * sind(angleofattack) for i in 1:N]

    # Get Pressure at NODES
    cpi = 1.0 .- vti .^ 2

    # Get Mean Pressure at PANEL MIDPOINTS
    cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0

    # Get Panel Vectors
    dn = [
        meshes[j].airfoil_nodes[2:end] .- meshes[j].airfoil_nodes[1:(end - 1)] for j in 1:M
    ][1]

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
    x0 = chord / 4.0
    z0 = 0.0

    #define portions of moment calculation
    cmmat = [2 1; 1 2] ./ 6.0
    dxddmi = [0.0 for i in 1:N]

    dxddmip1 = [0.0 for i in 1:N]

    offset = 0
    for m in 1:M
        for i in 1:(Ns[m] - 1)
            dxddmi[i + offset] =
                dn[i + offset][1] * (meshes[m].airfoil_nodes[i][1] - x0) +
                dn[i + offset][2] * (meshes[m].airfoil_nodes[i][2] - z0)

            dxddmip1[i + offset] =
                dn[i + offset][1] * (meshes[m].airfoil_nodes[i + 1][1] - x0) +
                dn[i + offset][2] * (meshes[m].airfoil_nodes[i + 1][2] - z0)
        end
        offset += Ns[m]
    end

    #calculate moment coefficient about leading edge
    cmi =
        sum([[cpi[i] cpi[i + 1]] * cmmat * [dxddmi[i]; dxddmip1[i]] for i in 1:(N - 1)]) / chord^2

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
