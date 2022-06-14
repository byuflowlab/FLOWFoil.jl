#=
Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 16 May 2022

Change Log:
=#

"""
    inviscid_polar(inviscid_solution, angleofattack; cascade=false)

Generate Polar object for inviscid system at given angle of attack.

**Arguements:**
 - `inviscid_solution::InviscidSolution` : Inviscid Solution object
 - `angleofattack::Float` : Angle of attack, in degrees
"""
function inviscid_polar(inviscid_solution, angleofattack; cascade=false)
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

    # Get Lift, Drag, and Moment Coefficients
    # leading edge location for moment reference)
    # _, leadingedgeidx = findmin(getindex.(nodes, 1))
    # x0 = nodes[leadingedgeidx][1]
    # z0 = nodes[leadingedgeidx][2]
    #quarter chord location (moment reference location for inviscid case)
    x0 = chord / 4.0
    z0 = 0.0

    # initialize pieces of moment calculation
    cmmat = [2 1; 1 2] ./ 6.0
    dxddmi = [0.0 for i in 1:N]
    dxddmip1 = [0.0 for i in 1:N]

    # initialize panel lengths
    dn = [[0.0 0.0] for i in 1:(N - 1)]

    offset = 0
    for m in 1:M
        for i in 1:(Ns[m] - 1)
            dxddmi[i + offset] =
                dn[i + offset][1] * (meshes[m].airfoil_nodes[i][1] - x0) +
                dn[i + offset][2] * (meshes[m].airfoil_nodes[i][2] - z0)

            dxddmip1[i + offset] =
                dn[i + offset][1] * (meshes[m].airfoil_nodes[i + 1][1] - x0) +
                dn[i + offset][2] * (meshes[m].airfoil_nodes[i + 1][2] - z0)
            dn[i + offset] = meshes[m].airfoil_nodes[i + 1] .- meshes[m].airfoil_nodes[i]
        end
        offset += Ns[m]
    end

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
 - `inviscid_solution::InviscidSolution` : the inviscid solution from which to find the vortex magnitudes.
 - `angleofattack::Float` : the angle of attack in degrees.
"""
function get_vortex_magnitudes(inviscid_solution, angleofattack)

    #rename for convenience
    gammas = inviscid_solution.panelgammas

    return [
        gammas[i, 1] * cosd(angleofattack) + gammas[i, 2] * sind(angleofattack) for
        i in 1:length(gammas)
    ]
end

"""
    function calculate_stream_grid(problem, solution, xrange, zrange; Nx=100, Nz=100)

Calculate stream function values across x- and z-ranges.

**Arguments:**
 - `problem::Problem` : Problem object associated with solution object.
 - `solution::Solution` : Solution object associated with problem.
 - `xrange::Array{Float}` : minimum and maximum values of grid in x-direction, [xmin xmax]
 - `zrange::Array{Float}` : minimum and maximum values of grid in z-direction, [zmin zmax]

**Keyword Arguments:**
 - `Nx::Int` : Number of gridpoints in x-direction
 - `Nz::Int` : Number of gridpoints in z-direction

**Returns:**
 - `xgrid::Array{Float}` : Array of x-locations of gridpoints
 - `zgrid::Array{Float}` : Array of z-locations of gridpoints
 - `streamgrid::Array{Float,2}` : Matrix of stream function values at gridpoints

"""
function calculate_stream_grid(problem, solution, xrange, zrange; Nx=100, Nz=100)

    #TODO: figure out dimensions throughout...
    vinf = 1.0

    # unpack for convenience
    meshes = solution.meshes
    nmeshes = length(meshes)
    N, Ns = size_system(meshes)
    offset = get_offset(Ns)
    alpha = problem.angleofattack
    panelgammas = solution.panelgammas
    nodes = solution.meshes[1].airfoil_nodes

    # get grid points
    xgrid = range(xrange[1]; stop=xrange[2], length=Nx)
    zgrid = range(zrange[1]; stop=zrange[2], length=Nz)

    # Calculate gamma on each panel
    gammas = get_gamma_magnitudes(panelgammas, alpha)

    # initialize grid
    streamgrid = [0.0 for z in zgrid, x in xgrid]

    # add in contributions for each airfoil
    for i in 1:nmeshes
        streamgrid .+= [
            get_stream_grid_value(;
                gammas=gammas[(1 + offset[i]):(Ns[i] + offset[i])],
                sigmas=nothing, #TODO: update this when viscous solver is done
                nodes=meshes[i].airfoil_nodes,
                point=[x z],
                vinf=vinf,
                alpha=alpha,
                blunt_te=meshes[i].blunt_te,
                txp=meshes[i].txp,
                tdp=meshes[i].tdp,
            ) for z in zgrid, x in xgrid
        ]
    end

    return xgrid, zgrid, streamgrid
end

"""
    get_stream_grid_value(gammas, nodes, point, vinf, alpha, blunt_te=false, txp=0.0, tdp=0.0)

For a given point, calculate the total stream function value induced by an airfoil.

**Arguments:**
 - `gammas::Array{Float}` : Array of vortex strength magnitudes at each of the nodes.
 - `sigmas::Array{Float}` : Array of source strength magnitudes at each of the nodes.
 - `nodes::Array{Array{Float}}` : Array of node locations
 - `point::Array{Float}` : [x z] location of the point in question
 - `vinf::Float` : Freestream magnitude
 - `alpah::Float` : Angle of attack in degrees
 - `blunt_te::Bool` : Flag whether airfoil has a blunt trailing edge or not
 - `txp::Float` : txp field from mesh object where nodes came from
 - `tdp::Float` : tdp field from mesh object where nodes came from
"""
function get_stream_grid_value(;
    gammas=nothing,
    sigmas=nothing,
    nodes=nothing,
    point=nothing,
    vinf=1.0,
    alpha=0.0,
    blunt_te=false,
    txp=0.0,
    tdp=0.0,
)

    # Start with freestream
    psi = vinf * (cosd(alpha) * point[2] - sind(alpha) * point[1])

    # Add contributions from each of the airfoil nodes
    for i in 1:(length(gammas) - 1)
        aij, aijp1 = get_vortex_influence(nodes[i], nodes[i + 1], point)
        psi += aij * gammas[i] + aijp1 * gammas[i + 1]
    end

    # add in trailing edge contributions
    if blunt_te
        psi +=
            0.5 *
            (gammas[end] - gammas[1]) *
            (
                txp * get_source_influence(nodes[end], nodes[1], point) -
                tdp * sum(get_vortex_influence(nodes[end], nodes[1], point))
            )
    end

    if sigmas != nothing
        @warn("Viscous solver not yet implemented")
        # Add in here the node source contributions
    end

    return psi
end

"""
    get_gamma_magnitudes(panelgammas, angleofattack)

Calculate vortex strength magnitudes for each of the nodes.

**Arguements:**
 - `panelgammas::Array{Float,2}` : panelgammas field from solution object.
 - `angleofattack::Float` : Angle of attack in degrees.
"""
function get_gamma_magnitudes(panelgammas, angleofattack)
    return [
        panelgammas[i, 1] * cosd(angleofattack) + panelgammas[i, 2] * sind(angleofattack)
        for i in 1:length(panelgammas[:, 1])
    ]
end
