# - If single airfoil, need to put Panel object in a vector - #
function generate_system_geometry(method::Xfoil, panel_geometry::NamedTuple)
    return generate_system_geometry(method, [panel_geometry])
end

"""
    generate_system_geometry(method::Xfoil, panel_geometry::AbstractVector; gap_tolerance=1e-10)

Constructs the system geometry for a multi-airfoil panel method simulation.
This includes indexing arrays, combined node and panel data, chord lengths, and precomputed influence geometry for each panel and trailing edge panel.

# Arguments
- `method::Xfoil`: The aerodynamic method object (currently unused, but kept for API consistency).
- `panel_geometry::AbstractVector`: Vector of panel geometry named tuples, one per airfoil/body, containing nodes, panel edges, vectors, lengths, etc.

# Keyword Arguments
- `gap_tolerance`: Threshold tolerance to detect blunt (non-sharp) trailing edges, default `1e-10`.

# Returns
- A named tuple `system_geometry` containing:
  - `nbodies`: Number of bodies/airfoils.
  - `panel_indices`: Vector of panel index ranges for each body.
  - `node_indices`: Vector of node index ranges for each body.
  - `nodes`: Combined matrix of all nodes for all bodies.
  - `mesh2panel`: Mapping from global mesh index to panel index.
  - `chord_length`: Overall chord length of the system (max trailing edge - min leading edge).
  - `panel_length`: Vector of panel lengths for all panels combined.
  - `r1`, `lnr1`, `r1normal`, `r1tangent`, `theta1`, `r2`, `lnr2`, `theta2`: Matrices of precomputed influence geometry quantities for all panels on all field points.
  - `TE_geometry`: Named tuple with trailing edge gap panel geometry and related influence parameters:
    - `blunt_te`: Boolean vector indicating blunt trailing edges.
    - `panel_length`: Trailing edge gap lengths.
    - `tdp`, `txp`: Trailing edge panel parameters.
    - `r1`, `lnr1`, `r1normal`, `r1tangent`, `theta1`, `r2`, `lnr2`, `theta2`: Influence geometry matrices for trailing edge panels.
"""
function generate_system_geometry(
    method::Xfoil, panel_geometry::AbstractVector; gap_tolerance=1e-10
)

    ### --- Convenience Variables --- ###
    nbodies = length(panel_geometry)
    npanels = [panel_geometry[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)
    nnodes = npanels .+ 1
    total_nodes = sum(nnodes)
    nodes = reduce(vcat, panel_geometry[i].nodes for i in 1:nbodies)

    # - Define Body Indexing - #

    #find starting indices for each body
    csnodes = cumsum(nnodes)
    cspanels = cumsum(npanels)

    node_indices = [(1 + (i == 1 ? 0 : csnodes[i - 1])):(csnodes[i]) for i in 1:nbodies]
    panel_indices = [(1 + (i == 1 ? 0 : cspanels[i - 1])):(cspanels[i]) for i in 1:nbodies]

    # - Map indices - #
    mesh2panel = reduce(vcat, [1:npanels[i] for i in 1:nbodies])

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([panel_geometry[i].panel_lengths[1] for i in 1:nbodies]))

    # Individual Chord Lengths
    leading_edges = zeros(TF, nbodies)
    trailing_edges = zeros(TF, nbodies)

    ### --- General Mesh Fields --- ###
    # Panel Length (contained in panels objects)
    panel_length = zeros(TF, (total_panels))

    # Distance from influencing panel edge 1 to field points (panel edges)
    r1 = zeros(TF, (total_nodes, total_panels))
    # # Distance from influencing panel edge 2 to field points (panel edges)
    r2 = zeros(TF, (total_nodes, total_panels))

    # Portion of distance 1 normal to the panel of influence
    r1normal = zeros(TF, (total_nodes, total_panels))
    # Portion of distance 1 tangent to the panel of influence
    r1tangent = zeros(TF, (total_nodes, total_panels))

    # Natural log of distance 1
    lnr1 = zeros(TF, (total_nodes, total_panels))
    # # Natural log of distance 2
    lnr2 = zeros(TF, (total_nodes, total_panels))

    # Angle between panel and field point with vertex at edge 1 of influencing panel
    theta1 = zeros(TF, (total_nodes, total_panels))
    # # Angle between panel and field point with vertex at edge 2 of influencing panel
    theta2 = zeros(TF, (total_nodes, total_panels))

    ### --- Trailing Edge Mesh Fields --- ###

    # Panel Length (contained in panels objects)
    TE_panel_length = zeros(TF, nbodies)

    # Distance from influencing panel edge 1 to field points (panel edges)
    r1_TE = zeros(TF, (total_nodes, nbodies))
    # # Distance from influencing panel edge 2 to field points (panel edges)
    r2_TE = zeros(TF, (total_nodes, nbodies))

    # Portion of distance 1 normal to the panel of influence
    r1normal_TE = zeros(TF, (total_nodes, nbodies))
    # Portion of distance 1 tangent to the panel of influence
    r1tangent_TE = zeros(TF, (total_nodes, nbodies))

    # Natural log of distance 1
    lnr1_TE = zeros(TF, (total_nodes, nbodies))
    # # Natural log of distance 2
    lnr2_TE = zeros(TF, (total_nodes, nbodies))

    # Angle between panel and field point with vertex at edge 1 of influencing panel
    theta1_TE = zeros(TF, (total_nodes, nbodies))
    # # Angle between panel and field point with vertex at edge 2 of influencing panel
    theta2_TE = zeros(TF, (total_nodes, nbodies))

    # Blunt Trailing Edge Flags
    blunt_te = zeros(Bool, nbodies)

    # tdp values
    tdp = zeros(TF, nbodies)

    # txp values
    txp = zeros(TF, nbodies)

    # trailing edge gap values
    trailing_edge_gap = zeros(TF, nbodies)

    ##### ----- Loop through each of the bodies to be influenced ----- #####
    for m in 1:nbodies

        ### --- Get Body Information --- ###

        # Get trailing edge information for each body
        tdp[m], txp[m], trailing_edge_gap[m], TE_panel_edges, TE_panel_vector, TE_panel_length[m] = get_trailing_edge_info(
            panel_geometry[m].panel_edges
        )

        # Get chord length
        leading_edges[m] = minimum(panel_geometry[m].panel_edges[:, :, 1])
        trailing_edges[m] = maximum(panel_geometry[m].panel_edges[:, :, 1])

        # Check if trailing edge is sharp or not
        if abs(TE_panel_length[m]) > gap_tolerance * (trailing_edges[m] - leading_edges[m])

            # set blunt_te to true
            blunt_te[m] = true

        else #(closed trailing edge)

            #set blunt_te to false
            blunt_te[m] = false
        end

        # Copy over panel lengths for convenience
        for pi in panel_indices
            panel_length[pi] = panel_geometry[m].panel_lengths
        end

        ##### ----- Loop through each of the bodies doing the influencing ----- #####
        for n in 1:nbodies

            ### --- Loop through each field node (panel edge) in body m --- ###
            for i in node_indices[m]

                # # Panel edges of panels being influenced (body m)
                # # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice
                # panidx = i > node_indices[m][end] -m ? i - m : i
                # field_panel_edge = panels[m].panel_edges[mesh2panel[panidx], :, :]

                # # Get vector and magnitude from first edge of the panel of influence to the field point (edge of panel being influenced)
                # # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice, using the second panel edge on the repeat
                # edgeidx = i == node_indices[m][end] ? 2 : 1

                ### --- Assemble Trailing Edge Gap Panel Infuences --- ###

                # Get trailing edge panel geometry

                # Calculate Influence Geometry
                r1_TE[i, m], r2_TE[i, m], r1normal_TE[i, m], r1tangent_TE[i, m], theta1_TE[i, m], theta2_TE[i, m], lnr1_TE[i, m], lnr2_TE[i, m] = calculate_influence_geometry(
                    TE_panel_edges,
                    TE_panel_vector,
                    trailing_edge_gap[m],
                    # field_panel_edge[edgeidx, :];
                    nodes[i, :];
                    gap_tolerance=gap_tolerance,
                )

                ### --- Loop through each influence panel in body n --- ###
                for j in panel_indices[n]

                    # - Rename For Convenience - #
                    # Panel edges of influencing panels (body n)
                    influence_panel_edge = panel_geometry[n].panel_edges[
                        mesh2panel[j], :, :,
                    ]

                    # Panel Vector and Length
                    influence_panel_vector = panel_geometry[n].panel_vectors[
                        mesh2panel[j], :,
                    ]

                    # - Calculate Influence Geometry - #

                    r1[i, j], r2[i, j], r1normal[i, j], r1tangent[i, j], theta1[i, j], theta2[i, j], lnr1[i, j], lnr2[i, j] = calculate_influence_geometry(
                        influence_panel_edge,
                        influence_panel_vector,
                        panel_length[mesh2panel[j]],
                        # field_panel_edge[edgeidx, :];
                        nodes[i, :];
                        gap_tolerance=gap_tolerance,
                    )
                end #for panels being influenced
            end #for influencing panels
        end #for body n
    end #for body m

    # Define System Chord Length
    # TODO: figure out how to expose different options here to the user in a slick way
    chord_length = maximum(trailing_edges) - minimum(leading_edges)

    # - Generate mesh objects - #
    #= TODO:
        There has got to be a better way to do all this.
        The problem is that it's unknown beforehand if there is a trailing edge gap.
        Possibly do the trailing edge gap panel discovery and geometry in the paneling function(s), and then pass those in as part of everything.
        Indexing could get tricky trying to combine everything together though, since the TE panels only influence and are not influenced, so they don't add any more equations to the system.
    =#

    # trailing edge gap panel geomtry
    TE_geometry = (;
        blunt_te,
        panel_length=trailing_edge_gap,
        tdp,
        txp,
        r1=r1_TE,
        lnr1=lnr1_TE,
        r1normal=r1normal_TE,
        r1tangent=r1tangent_TE,
        theta1=theta1_TE,
        r2=r2_TE,
        lnr2=lnr2_TE,
        theta2=theta2_TE,
    )

    # main system_geometry
    system_geometry = (;
        nbodies,
        panel_indices,
        node_indices,
        nodes,
        mesh2panel,
        chord_length,
        panel_length,
        r1,
        lnr1,
        r1normal,
        r1tangent,
        theta1,
        r2,
        lnr2,
        theta2,
        TE_geometry,
    )

    return system_geometry
end

"""
    calculate_influence_geometry(
        influence_panel_edge,
        influence_panel_vector,
        influence_panel_length,
        field_panel_edge;
        gap_tolerance=1e-10,
    )

Calculate all the various geometric pieces of the influence coefficients.

# Arguments
- `influence_panel_edge::Array{Float}` : Array of the edge locations of the panel doing the influencing.
- `influence_panel_vector::Vector{Float}` : Vector from the first to second panel edge.
- `influence_panel_length::Float` : Length of the influencing panel.
- `field_point::Vector::Float` : position of field point (point being influence).

# Keyword Arguments
- `gap_tolerance::Float` : Distance at which self-induction is assumed (to avoid NaNs and Infs). default = 1e-10

# Returns
- `r1::Matrix{TF}` : distance from first panel edge of influence panel to field point
- `r2::Matrix{TF}` : distance from second panel edge of influence panel to field point.
- `r1normal::Matrix{TF}` : component of r1 in normal direction of influence panel
- `r1tangent::Matrix{TF}` : component of r1 in tangent direction of influence panel
- `theta1::Matrix{TF}` : angle between influece panel vector and the r1 vector, adjusted to π for self-influence of panels
- `theta2::Matrix{TF}` angle between influence panel vector and the r2 vector, adjusted to π for self-influence of panels
- `lnr1::Matrix{TF}` : ln(r1), adjusted to zero for self-influence of panels
- `lnr2::Matrix{TF}` : ln(r2), adjusted to zero for self-influence of panels
"""
function calculate_influence_geometry(
    influence_panel_edge,
    influence_panel_vector,
    influence_panel_length,
    field_point;
    gap_tolerance=1e-10,
)

    # Get vector and magnitude from first edge of the panel of influence to the field point (edge of panel being influenced)
    r1vec, r1 = get_r(influence_panel_edge[1, :], field_point)

    # Get vector and magnitude from second edge of the panel of influence to the field point (edge of panel being influenced)
    # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice, using the second panel edge on the repeat
    r2vec, r2 = get_r(influence_panel_edge[2, :], field_point)

    # Get the component of r1vec normal to the panel of influence
    r1normal = get_r_normal(r1vec, influence_panel_vector, influence_panel_length)

    # Get the component of r1vec tangent to the panel of influence
    r1tangent = get_r_tangent(r1vec, influence_panel_vector, influence_panel_length)

    # Get the angle between the panel of influence and field point with angle vertex at the first edge of the panel of influence
    theta1 = get_theta(r1normal, r1tangent)

    # Get the angle between the panel of influence and field point with angle vertex at the second edge of the panel of influence
    theta2 = get_theta2(r1normal, r1tangent, influence_panel_length)

    # Calculate natural log values
    #= Natural Log values are to be used in calculating the influence coefficients
        We set the natural logs to zero if the field point is very close to the panel,
        that is, for the self-influence cases.
        NOTE: probably will have other probems if another panel is as close as a panel is to itself...
        We also adjust the values of the angles previously calculated for self-induction cases to be pi or zero depending on positions.
    =#
    if r1 < gap_tolerance
        lnr1 = 0.0
        theta1 = pi
        theta2 = pi

    else
        lnr1 = log(r1)
    end #if

    if r2 < gap_tolerance
        lnr2 = 0.0
        theta1 = 0.0
        theta2 = 0.0

    else
        lnr2 = log(r2)
    end #if

    return r1, r2, r1normal, r1tangent, theta1, theta2, lnr1, lnr2
end
