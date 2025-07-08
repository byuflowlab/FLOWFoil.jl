"""
    get_gamma_magnitudes(strengths, angle_of_attack)

Calculate vortex strength magnitudes for each of the nodes.

# Arguments
 - `strengths::Array{Float,2}` : strengths field from solution object.
 - `angle_of_attack::Float` : Angle of attack in degrees.
"""
function get_gamma_magnitudes(strengths, angle_of_attack)
    return [
        strengths[i, 1] * cosd(angle_of_attack) + strengths[i, 2] * sind(angle_of_attack)
        for i in 1:length(strengths[:, 1])
    ]
end

function get_vortex_influence(
    panel_length, r1, r2, r1normal, r1tangent, theta1, theta2, lnr1, lnr2
)

    # get psibargamma value
    psibargamma = get_psibargamma(
        theta1, theta2, lnr1, lnr2, panel_length, r1normal, r1tangent
    )

    # get psitildegamma value
    psitildegamma = get_psitildegamma(
        psibargamma, r1, r2, theta1, theta2, lnr1, lnr2, panel_length, r1normal, r1tangent
    )

    # put psi`s together
    return (psibargamma - psitildegamma), psitildegamma
end

"""
    get_stream_grid_value(strengths, nodes, field_point, vinf, angle_of_attack, blunt_te=false, txp=0.0, tdp=0.0)

For a given field_point, calculate the total stream function value induced by an airfoil.

**Arguments:**
 - `strengths::Array{Float}` : Array of vortex strength magnitudes at each of the nodes.
 - `sigmas::Array{Float}` : Array of source strength magnitudes at each of the nodes.
 - `nodes::Array{Array{Float}}` : Array of node locations
 - `field_point::Array{Float}` : [x y] location of the field_point in question
 - `vinf::Float` : Freestream magnitude
 - `alpah::Float` : Angle of attack in degrees
 - `blunt_te::Bool` : Flag whether airfoil has a blunt trailing edge or not
 - `txp::Float` : txp field from mesh object where nodes came from
 - `tdp::Float` : tdp field from mesh object where nodes came from
"""
function get_stream_grid_value(
    panel_geometry,
    strengths,
    field_point;
    sigmas=nothing,
    vinf=1.0,
    angle_of_attack=0.0,
    blunt_te=false,
    txp=0.0,
    tdp=0.0,
)

    # Start with freestream
    psi =
        vinf *
        (cosd(angle_of_attack) * field_point[2] - sind(angle_of_attack) * field_point[1])

    # Add contributions from each of the airfoil nodes
    for i in 1:(length(strengths) - 1)
        aij, aijp1 = get_vortex_influence(
            panel_geometry.panel_lengths[i],
            calculate_influence_geometry(
                panel_geometry.panel_edges[i, :, :],
                panel_geometry.panel_vectors[i, :],
                panel_geometry.panel_lengths[i],
                field_point,
            )...,
        )

        psi += aij * strengths[i] + aijp1 * strengths[i + 1]
    end

    # add in trailing edge contributions
    if blunt_te
        psi +=
            0.5 *
            (strengths[end] - strengths[1]) *
            (
                txp * get_source_influence(nodes[end], nodes[1], field_point) -
                tdp * sum(get_vortex_influence(nodes[end], nodes[1], field_point))
            )
    end

    if sigmas != nothing
        @warn("Viscous solver not yet implemented")
        # Add in here the node source contributions
    end

    return psi
end

"""
    function calculate_stream_grid(problem, solution, xrange, zrange; Nx=100, Nz=100)

Calculate stream function values across x- and y-ranges.

# Arguments
 - `xrange::Array{Float}` : minimum and maximum values of grid in x-direction, [xmin xmax]
 - `zrange::Array{Float}` : minimum and maximum values of grid in y-direction, [zmin zmax]

# Keyword Arguments
 - `Nx::Int` : Number of gridpoints in x-direction
 - `Nz::Int` : Number of gridpoints in y-direction

# Returns
 - `xgrid::Array{Float}` : Array of x-locations of gridpoints
 - `zgrid::Array{Float}` : Array of y-locations of gridpoints
 - `stream_grid::Array{Float,2}` : Matrix of stream function values at gridpoints

"""
function calculate_stream_grid(
    method::Xfoil, nbodies, panel_geometry, system_geometry, strengths
)

    # unpack for convenience
    nodes = panel_geometry.nodes

    # get grid points
    xgrid = range(method.xrange[1]; stop=method.xrange[2], length=method.Nx)
    zgrid = range(method.zrange[1]; stop=method.zrange[2], length=method.Nz)

    # Calculate gamma on each panel
    strengths = get_gamma_magnitudes(strengths, method.angle_of_attack)

    # initialize grid
    stream_grid = [0.0 for y in zgrid, x in xgrid]

    # add in contributions for each airfoil
    nidx = system_geometry.node_indices
    for m in 1:nbodies
        stream_grid .+= [
            get_stream_grid_value(
                panel_geometry,
                strengths[nidx[m], :],
                [x; y];
                sigmas=nothing, #TODO: update this when viscous solver is done
                vinf=1.0,
                angle_of_attack=method.angle_of_attack,
                blunt_te=system_geometry.TE_geometry.blunt_te[m],
                txp=system_geometry.TE_geometry.txp[m],
                tdp=system_geometry.TE_geometry.tdp[m],
            ) for y in zgrid, x in xgrid
        ]
    end

    return xgrid, zgrid, stream_grid
end
