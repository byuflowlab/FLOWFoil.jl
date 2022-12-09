#=
Planar Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 16 May 2022

Change Log:
=#

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
    nodes = solution.meshes[1].nodes

    # get grid points
    xgrid = range(xrange[1]; stop=xrange[2], length=Nx)
    zgrid = range(zrange[1]; stop=zrange[2], length=Nz)

    # Calculate gamma on each panel
    gammas = get_gamma_magnitudes(panelgammas, alpha)

    # initialize grid with freestream
    streamgrid = [vinf * (cosd(alpha) * z - sind(alpha) * x) for z in zgrid, x in xgrid]

    # add in contributions for each airfoil
    for i in 1:nmeshes
        streamgrid .+= [
            get_stream_grid_value(;
                gammas=gammas[(1 + offset[i]):(Ns[i] + offset[i])],
                sigmas=nothing, #TODO: update this when viscous solver is done
                nodes=meshes[i].nodes,
                point=[x z],
                vinf=vinf,
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
    blunt_te=false,
    txp=0.0,
    tdp=0.0,
)

    #initialize
    psi = 0.0

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
