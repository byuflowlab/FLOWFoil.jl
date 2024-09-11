"""
    MfoilOutputs{TF} <: Outputs

**Fields:**
- `lift::Matrix{Float}` : Lift Coefficient.
- `drag::Matrix{Float}` : Total Drag Coefficient.
- `pdrag::Matrix{Float}` : Pressure Drag Coefficient.
- `idrag::Matrix{Float}` : Induced Drag Coefficient.
- `moment::Matrix{Float}` : Moment Coefficient.
- `surfacevelocity::Array{Float}` : smoothed surface velocity distribution
- `surfacepressure::Array{Float}` : smoothed surface pressure distribution
- `xsmooth::Array{Float}` : x-values associated with smoothed surface distributions
"""
struct MfoilOutputs{TF} <: Outputs
    lift::Matrix{TF}
    drag::Matrix{TF}
    pdrag::Matrix{TF}
    idrag::Matrix{TF}
    moment::Matrix{TF}
    surface_velocity::Array{TF,3}
    surface_pressure::Array{TF,3}
    xsmooth::Array{TF,3}
end

function post_process(
    ::Mfoil, problem, panels, mesh, solution; npanels=80, debug=false
)

    ##### ----- Set Up ----- #####

    ### --- Rename for Convenience --- ###

    # number of airfoils
    nbodies = mesh.nbodies

    # number of angles of attack
    flow_angle = problem.flow_angle
    naoa = length(flow_angle)

    # node indices
    pidx = mesh.panel_indices
    nidx = mesh.node_indices

    # vortex strengths
    gamma0 = solution.x[:, 1]
    gamma90 = solution.x[:, 2]

    # chord length
    chord = mesh.chord

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
    if debug
        # Surface Velocities
        v_surf = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
        # Surface Pressures
        p_surf = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
    else
        # Surface Velocities
        v_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
        # Surface Pressures
        p_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    end

    smooth_nodes = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    ##### ----- Loop Through Bodies ----- #####
    for m in 1:nbodies

        # Panel Values
        # length
        panel_length = panels[m].panel_length
        # vector
        panel_vector = panels[m].panel_vector
        # edge locations
        panel_edges = panels[m].panel_edges
        # midpoints
        panel_midpoints = (panel_edges[:, 2, :] .+ panel_edges[:, 1, :]) ./ panel_length

        for a in 1:naoa

            ### --- Get Surface Distributions --- ###
            # - Get Surface Velocity - #
            #= NOTE:
                For xfoil-like method, the "solution" values ARE the velocity components on the surface.
            =#
            vti = [
                gamma0[i] * cosd(flow_angle[a]) + gamma90[i] * sind(flow_angle[a]) for
                i in nidx[m]
            ]

            # - Get Surface Pressure (Steady State) - #
            cpi = 1.0 .- vti .^ 2

            # Get Mean Pressure at PANEL MIDPOINTS
            cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0

            # organize outputs based on debug flag
            if debug
                v_surf[m][:, a] = vti
                p_surf[m][:, a] = cpi
                # Get Mean Pressure at PANEL MIDPOINTS
                cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0
            else
                #smooth_distributions functions are found in utils.jl
                v_surf[m, :, a], smooth_nodes[m, :, a] = smooth_distributions(
                    Linear(), panel_edges, vti, npanels
                )
                p_surf[m, :, a], _ = smooth_distributions(
                    Linear(), panel_edges, cpi, npanels
                )
            end

            #---------------------------------#
            #      Calculate Coefficients     #
            #---------------------------------#

            #quarter chord location (moment reference location for inviscid case)
            x0 = chord / 4.0
            z0 = 0.0 #x0*sind(flow_angle[a]) #TODO should this be zero, or rotated with the airfoil?

            panelidx = mesh.mesh2panel

            ### --- Calculate Lift Coefficient --- ###
            cl[m, a] =
                sum([
                    cpibar[i] * (
                        -sind(flow_angle[a]) * panel_vector[i, 2] -
                        cosd(flow_angle[a]) * panel_vector[i, 1]
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
                        cosd(flow_angle[a]) * panel_vector[i, 2] -
                        sind(flow_angle[a]) * panel_vector[i, 1]
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

    # Again, returns depend on debug flag
    if debug
        return cl, cd, cdp, cdi, cm, v_surf, p_surf
    else
        return MfoilOutputs(cl, cd, cdp, cdi, cm, v_surf, p_surf, smooth_nodes)
    end
end
