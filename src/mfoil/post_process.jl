"""
    post_process(method::Mfoil, panel_geometry::NamedTuple, system_geometry, strengths, flow_angles; npanels=80)

Post-processes the aerodynamic solution for a single airfoil by wrapping the vector version of `post_process`.

# Arguments
- `method::Mfoil`: The Mfoil method object (configuration).
- `panel_geometry::NamedTuple`: Panel geometry data for a single airfoil.
- `system_geometry`: System geometry containing info like panel indices and chord lengths.
- `strengths`: Matrix of vortex strengths (rows correspond to nodes/panels, columns to different components).
- `flow_angles`: Vector of flow angles of attack (in degrees).

# Keyword Arguments
- `npanels`: Number of panels to use for output (default 80).

# Returns
- An `InviscidOutputs` object containing aerodynamic coefficients and surface distributions.
"""
function post_process(
    method::Mfoil,
    panel_geometry::NamedTuple,
    system_geometry,
    strengths,
    flow_angles;
    npanels=80,
)
    return post_process(method, [panel_geometry], system_geometry, strengths, flow_angles)
end

"""
    post_process(::Mfoil, panel_geometry::AbstractVector, system_geometry, strengths, flow_angles)

Computes aerodynamic surface velocities, pressures, and coefficients (lift, drag, moment) for one or more airfoils over a range of flow angles.

# Arguments
- `::Mfoil`: Method configuration object (not used directly but required for dispatch).
- `panel_geometry::AbstractVector`: Vector of panel geometry data, one for each airfoil/body.
- `system_geometry`: System geometry object containing details such as panel and node indices, chord lengths, and mapping arrays.
- `strengths`: Matrix containing vortex strengths at nodes for each angle of attack.
- `flow_angles`: Vector of flow angles of attack in degrees.

# Returns
- If one airfoil, returns an `InviscidOutputs` struct with fields:
  - `vs`: Surface velocity distributions.
  - `cp`: Surface pressure distributions.
  - `cl`: Lift coefficient array (angles × bodies).
  - `cd`: Total drag coefficient array.
  - `cm`: Moment coefficient array.
- If multiple airfoils, returns an `InviscidOutputs` struct with surface velocities and pressures as vectors corresponding to each body.
"""
function post_process(
    ::Mfoil,
    panel_geometry::AbstractVector,
    system_geometry,
    strengths,
    flow_angles; #debug=false
)

    ##### ----- Set Up ----- #####

    ### --- Rename for Convenience --- ###

    # number of airfoils
    nbodies = system_geometry.nbodies

    # number of angles of attack
    naoa = length(flow_angles)

    # node indices
    pidx = system_geometry.panel_indices
    nidx = system_geometry.node_indices

    # vortex strengths
    gamma0 = strengths[:, 1]
    gamma90 = strengths[:, 2]

    # chord length
    chord = system_geometry.chord_length

    ### --- Initialize Outputs --- ###

    # output floating point type
    TF = typeof(chord)

    # - Coefficients - #

    # Lift coefficient
    cl = zeros(TF, naoa, nbodies)

    # Total drag coefficient
    cd = zeros(TF, naoa, nbodies)
    # Profile drag coefficient
    cdp = zeros(TF, naoa, nbodies)
    # Inviscid drag coefficient
    cdi = zeros(TF, naoa, nbodies)

    # Moment coefficient
    cm = zeros(TF, naoa, nbodies)

    # Surface Velocities
    vs = [zeros(nidx[m][end] - nidx[m][1] + 1, naoa) for m in 1:nbodies]

    # Surface Pressures
    cp = [zeros(nidx[m][end] - nidx[m][1] + 1, naoa) for m in 1:nbodies]

    ##### ----- Loop Through Bodies ----- #####
    for m in 1:nbodies

        # Panel Values
        # length
        # panel_lengths = panel_geometry[m].panel_lengths
        # vector
        panel_vector = panel_geometry[m].panel_vectors
        # edge locations
        panel_edges = panel_geometry[m].panel_edges
        # midpoints
        # panel_midpoints = (panel_edges[:, 2, :] .+ panel_edges[:, 1, :]) ./ panel_lengths

        for a in 1:naoa

            ### --- Get Surface Distributions --- ###
            # - Get Surface Velocity - #
            #= NOTE:
                For xfoil-like method, the "solution" values ARE the velocity components on the surface.
            =#
            vti = [
                gamma0[i] * cosd(flow_angles[a]) + gamma90[i] * sind(flow_angles[a]) for
                i in nidx[m]
            ]

            # - Get Surface Pressure (Steady State) - #
            cpi = 1.0 .- vti .^ 2

            # Get Mean Pressure at PANEL MIDPOINTS
            cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0

            # - Organize surface velocities and pressures - #
            vs[m][:, a] = vti
            cp[m][:, a] = cpi

            #---------------------------------#
            #      Calculate Coefficients     #
            #---------------------------------#

            #quarter chord location (moment reference location for inviscid case)
            x0 = chord / 4.0
            z0 = 0.0 #x0*sind(flow_angles[a]) #TODO should this be zero, or rotated with the airfoil?

            panelidx = system_geometry.mesh2panel

            ### --- Calculate Lift Coefficient --- ###
            cl[a, m] =
                sum([
                    cpibar[i] * (
                        -sind(flow_angles[a]) * panel_vector[i, 2] -
                        cosd(flow_angles[a]) * panel_vector[i, 1]
                    ) for i in panelidx[pidx[m]]
                ]) / chord

            ### --- Calculate Drag Coefficients --- ###
            #= NOTE:
                For the inviscid case, cdp is zero and cdi=cd for now.
            =#
            cdp[a, m] = 0.0
            #TODO: check of drag calculation needs to be negative or not
            cdi[a, m] =
                -sum([
                    cpibar[i] * (
                        cosd(flow_angles[a]) * panel_vector[i, 2] -
                        sind(flow_angles[a]) * panel_vector[i, 1]
                    ) for i in panelidx[pidx[m]]
                ]) / chord
            cd[a, m] = cdi[a, m]

            ### --- Calculate Moment Coefficient --- ###
            # initialize pieces of moment calculation
            cmmat = [2 1; 1 2] ./ 6.0

            # Moment arms
            dxddmi =
                panel_vector[panelidx[pidx[m]], 1] .*
                (panel_edges[panelidx[pidx[m]], 1, 1] .- x0) .+
                panel_vector[panelidx[pidx[m]], 2] .*
                (panel_edges[panelidx[pidx[m]], 1, 2] .- z0)

            dxddmip1 =
                panel_vector[panelidx[pidx[m]], 1] .*
                (panel_edges[panelidx[pidx[m]], 2, 1] .- x0) .+
                panel_vector[panelidx[pidx[m]], 2] .*
                (panel_edges[panelidx[pidx[m]], 2, 2] .- z0)

            # Moment Coefficient Calculation
            cm[a, m] =
                sum([
                    ([cpi[i] cpi[i + 1]] * cmmat * [
                        dxddmi[i]
                        dxddmip1[i]
                    ])[1] for i in panelidx[pidx[m]]
                ]) / chord^2
        end
    end
    if nbodies == 1
        return InviscidOutputs(vs[1], cp[1], cl, cd, cm)
    else
        return InviscidOutputs(vs, cp, cl, cd, cm)
    end
end
