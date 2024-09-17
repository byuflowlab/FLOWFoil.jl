function post_process(method::Mfoil, panel_geometry::NamedTuple, system_geometry, strengths, flow_angles; npanels=80)
    return post_process(
        method, [panel_geometry], system_geometry, strengths, flow_angles
    )
end

function post_process(
    ::Mfoil, panel_geometry::AbstractVector, system_geometry, strengths, flow_angles; #debug=false
)

    ##### ----- Set Up ----- #####

    ### --- Rename for Convenience --- ###

    # number of airfoils
    nbodies = system_geometry.nbodies

    # number of panels
    npanels = panel_geometry[1].npanels

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
    cl = zeros(TF, nbodies, naoa)

    # Total drag coefficient
    cd = zeros(TF, nbodies, naoa)
    # Profile drag coefficient
    cdp = zeros(TF, nbodies, naoa)
    # Inviscid drag coefficient
    cdi = zeros(TF, nbodies, naoa)

    # Moment coefficient
    cm = zeros(TF, nbodies, naoa)

    # - Surface Distributions - #

    #= NOTE:
        If debug is true, then we can allow for the raw outputs for surface distributions which may or may not be the same length (thus requiring a tuple).
        If debug is false, then we assume that we'll be smoothing the outputs with Akima splines.
    =#
    # if debug
    #     # Surface Velocities
    #     tangential_velocities = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
    #     # Surface Pressures
    #     surface_pressures = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
    # else
        # Surface Velocities
        tangential_velocities = zeros(TF, nbodies, 2 * npanels - 1, naoa)
        # Surface Pressures
        surface_pressures = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    # end

    smooth_nodes = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    ##### ----- Loop Through Bodies ----- #####
    for m in 1:nbodies

        # Panel Values
        # length
        panel_lengths = panel_geometry[m].panel_lengths
        # vector
        panel_vector = panel_geometry[m].panel_vectors
        # edge locations
        panel_edges = panel_geometry[m].panel_edges
        # midpoints
        panel_midpoints = (panel_edges[:, 2, :] .+ panel_edges[:, 1, :]) ./ panel_lengths

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

            # organize outputs based on debug flag
            # if debug
            #     tangential_velocities[m][:, a] = vti
            #     surface_pressures[m][:, a] = cpi
            #     # Get Mean Pressure at PANEL MIDPOINTS
            #     cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0
            # else
                #smooth_distributions functions are found in utils.jl
                tangential_velocities[m, :, a], smooth_nodes[m, :, a] = smooth_distributions(
                    Mfoil(), panel_edges, vti, npanels
                )
                surface_pressures[m, :, a], _ = smooth_distributions(
                    Mfoil(), panel_edges, cpi, npanels
                )
                # Get Mean Pressure at PANEL MIDPOINTS
                cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0
            # end

            #---------------------------------#
            #      Calculate Coefficients     #
            #---------------------------------#

            #quarter chord location (moment reference location for inviscid case)
            x0 = chord / 4.0
            z0 = 0.0 #x0*sind(flow_angles[a]) #TODO should this be zero, or rotated with the airfoil?

            panelidx = system_geometry.mesh2panel

            ### --- Calculate Lift Coefficient --- ###
            cl[m, a] =
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
            cdp[m, a] = 0.0
            #TODO: check of drag calculation needs to be negative or not
            cdi[m, a] =
                -sum([
                    cpibar[i] * (
                        cosd(flow_angles[a]) * panel_vector[i, 2] -
                        sind(flow_angles[a]) * panel_vector[i, 1]
                    ) for i in panelidx[pidx[m]]
                ]) / chord
            cd[m, a] = cdi[m, a]

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
            cm[m, a] =
                sum([
                    ([cpi[i] cpi[i + 1]] * cmmat * [
                        dxddmi[i]
                        dxddmip1[i]
                    ])[1] for i in panelidx[pidx[m]]
                ]) / chord^2
        end
    end

    return (; cl, cd, cdp, cdi, cm, tangential_velocities, surface_pressures, smooth_nodes)
end
