#=
Planar Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 16 May 2022

Change Log:
=#

"""
    get_planar_polar(inviscid_solution, angleofattack; cascade=false)

Generate PlanarPolar object for inviscid system at given angle of attack.

**Arguements:**
 - `inviscid_solution::InviscidSolution` : Inviscid Solution object
 - `angleofattack::Float` : Angle of attack, in degrees
"""
function get_planar_polar(inviscid_solution, angleofattack::T;
                            chord=nothing, onlymeshes=nothing) where {T}
    M = length(inviscid_solution.Ns)            # Number of meshes

    # error case
    if onlymeshes != nothing
        for mi in onlymeshes
            @assert mi<=M "Requested mesh $mi, but max is $M"
        end
    end

    meshes = onlymeshes==nothing ?  inviscid_solution.meshes :
                                    inviscid_solution.meshes[onlymeshes]

    rtype = promote_type(T, eltype.(mesh.nodes[1] for mesh in meshes)...)

    # determine indices of panels in the meshes to be analyzed
    rngall = vcat(0, cumsum(inviscid_solution.Ns))
    rng = Iterators.flatten(( rngall[mi]+1:rngall[mi+1] for mi in (onlymeshes==nothing ? (1:M) : onlymeshes) ))

    # rename fields for convenience
    gamma0 = [inviscid_solution.panelgammas[i, 1] for i in rng]
    gamma90 = [inviscid_solution.panelgammas[i, 2] for i in rng]
    N = sum(1 for i in rng)

    # Minimum x position
    xmin = minimum(minimum(node[1] for node in mesh.nodes) for mesh in meshes)

    #calculate chord
    if chord != nothing
        c = chord
    else
        c = maximum(maximum(node[1] for node in mesh.nodes) for mesh in meshes) - xmin
    end

    # Get Velocity at NODES
    vti = gamma0*cosd(angleofattack) .+ gamma90*sind(angleofattack)

    # Get Pressure at NODES
    cpi = 1.0 .- vti.^2

    # Get Mean Pressure at PANEL MIDPOINTS
    cpibar = (cpi[1:(end-1)] .+ cpi[2:end]) ./ 2.0

    # Get Lift, Drag, and Moment Coefficients
    # leading edge location for moment reference)
    # _, leadingedgeidx = findmin(getindex.(nodes, 1))
    # x0 = nodes[leadingedgeidx][1]
    # z0 = nodes[leadingedgeidx][2]

    #quarter c location (moment reference location for inviscid case)
    x0 = xmin + c/4
    z0 = 0.0 #chord*sind(angleofattack) #TODO should this be zero, or rotated with the airfoil?

    # initialize pieces of moment calculation
    cmmat = [2 1; 1 2] ./ 6.0
    dxddmi = zeros(rtype, N)
    dxddmip1 = zeros(rtype, N)

    # initialize panel lengths
    dn = zeros(rtype, 2, N-1)

    p = 1                           # Panel count
    for mesh in meshes

        nn = length(mesh.nodes)
        for (node, nodep1) in zip(view(mesh.nodes, 1:nn-1), view(mesh.nodes, 2:nn))

            for i in 1:2
                dn[i, p] = nodep1[i]
                dn[i, p] -= node[i]
            end

            dxddmi[p] = dn[1, p]*(node[1] - x0) + dn[2, p]*(node[2] - z0)
            dxddmip1[p] = dn[1, p]*(nodep1[1] - x0) + dn[2, p]*(nodep1[2] - z0)

            p += 1
        end

    end

    # get lift coefficient
    sina, cosa = sind(angleofattack), cosd(angleofattack)
    cl = sum( cp*(-sina*d[2] - cosa*d[1]) for (cp, d) in zip(cpibar, eachcol(dn)) ) / c

    # get drag coefficient
    cdp = zero(rtype)
    cdi = sum( cp*(cosa*d[2] - sina*d[1]) for (cp, d) in zip(cpibar, eachcol(dn)) ) / c
    cd = cdi

    # calculate moment coefficient about leading edge
    ncp = length(cpi)
    cmi = sum( [cp cpp1] * cmmat * [dxddm; dxddmp1]
                        for (cp, cpp1, dxddm, dxddmp1)
                        in zip(view(cpi, 1:ncp-1), view(cpi, 2:ncp), dxddmi, dxddmip1)
            ) / c^2

    # Create PlanarPolar Object
    planar_polar = PlanarPolar(cl, cd, cdp, cdi, cmi[1], vti, cpi)

    return planar_polar
end

const inviscid_polar = get_planar_polar

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
