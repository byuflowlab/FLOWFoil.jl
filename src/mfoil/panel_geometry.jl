"""
"""function generate_panels(p::PlanarProblem, coordinates)

    #broadcast for multiple airfoils
    return reduce(vcat, generate_panels.(Ref(p), coordinates))
end

"""
"""
function generate_panels(::PlanarProblem, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    numpanels = length(x) - 1

    # Initialize Outputs
    panel_edges = zeros(TF, numpanels, 2, 2)
    panel_lengths = zeros(TF, numpanels)
    panel_vectors = zeros(TF, numpanels, 2)
    nodes = zeros(TF, numpanels + 1, 2)
    nodes[1, :] = [x[1] y[1]]
    nodes[end, :] = [x[end] y[end]]

    for i in 1:numpanels
        # Get node locations from x,y coordinates
        panel_edges[i, :, :] = [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panel_lengths[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panel_vectors[i, :] = [x[i + 1] - x[i] y[i + 1] - y[i]]

        # Set Node Values
        nodes[i, :] = [x[i] y[i]]
    end

    return [LinearFlatPanel(numpanels, panel_edges, panel_vectors, panel_lengths, nodes)]
end

"""
"""
function generate_panels!(::PlanarProblem, panels, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    numpanels = length(x) - 1

    # Initialize Outputs
    panels.nodes[1, :] .= [x[1] y[1]]
    panels.nodes[end, :] .= [x[end] y[end]]

    for i in 1:numpanels
        # Get node locations from x,y coordinates
        panels.panel_edges[i, :, :] .= [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panels.panel_lengths[i] .= sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panels.panel_vectors[i, :] .= [x[i + 1] - x[i] y[i + 1] - y[i]]

        # Set Node Values
        panels.nodes[i, :] .= [x[i] y[i]]
    end

    return nothing
end


