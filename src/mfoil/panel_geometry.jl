function generate_panel_geometry(method::Mfoil, coordinates)

    #broadcast for multiple airfoils
    return reduce(vcat, generate_panel_geometry.(Ref(method), coordinates))
end

function generate_panel_geometry(method::Mfoil, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    npanels = length(x) - 1

    # Initialize Outputs
    panel_geometry = (;
        # - Rename for Convenience - #
        npanels=npanels,
        # - Initialize Outputs - #
        panel_edges=zeros(TF, npanels, 2, 2),
        panel_vectors=zeros(TF, npanels, 2),
        panel_lengths=zeros(TF, npanels),
        nodes=zeros(TF, npanels + 1, 2)
    )

    return generate_panel_geometry!(method, panel_geometry, coordinates)
end

function generate_panel_geometry!(::Mfoil, panel_geometry, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Initialize Outputs
    panel_geometry.nodes[1, :] = [x[1]; y[1]]
    panel_geometry.nodes[end, :] = [x[end]; y[end]]

    for i in 1:panel_geometry.npanels
        # Get node locations from x,y coordinates
        panel_geometry.panel_edges[i, :, :] .= [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panel_geometry.panel_lengths[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panel_geometry.panel_vectors[i, :] = [x[i + 1] - x[i]; y[i + 1] - y[i]]

        # Set Node Values
        panel_geometry.nodes[i, :] = [x[i]; y[i]]
    end

    return panel_geometry
end
