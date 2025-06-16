"""
    generate_panel_geometry(method::Mfoil, coordinates)

Generates panel geometries for multiple airfoils by broadcasting over a collection of coordinate sets.

# Arguments
- `method::Mfoil`: The Mfoil method configuration object.
- `coordinates`: A collection (e.g. vector) of coordinate matrices, each representing airfoil node coordinates.

# Returns
- A concatenated vector of panel geometry objects for all provided airfoils.
"""
function generate_panel_geometry(method::Mfoil, coordinates)

    #broadcast for multiple airfoils
    return reduce(vcat, generate_panel_geometry.(Ref(method), coordinates))
end

"""
    generate_panel_geometry(method::Mfoil, coordinates::Matrix{TF}) where {TF}

Generates the panel geometry for a single airfoil given node coordinates.

# Arguments
- `method::Mfoil`: The Mfoil method configuration object.
- `coordinates::Matrix{TF}`: Nx2 matrix containing the x and y coordinates of the airfoil nodes, where N is the number of nodes.

# Returns
- Named tuple `panel_geometry` containing:
  - `npanels::Int`: Number of panels (nodes - 1).
  - `panel_edges::Array{TF, 3}`: An (npanels, 2, 2) array containing the start and end points of each panel.
  - `panel_vectors::Array{TF, 2}`: An (npanels, 2) array of vectors describing panel directions.
  - `panel_lengths::Array{TF}`: Vector of panel lengths.
  - `nodes::Array{TF, 2}`: Array of node coordinates (npanels + 1, 2).
"""
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

"""
    generate_panel_geometry!(method::Mfoil, panel_geometry, coordinates::Matrix{TF}) where {TF}

Fills the `panel_geometry` structure with geometric data calculated from airfoil node coordinates.

# Arguments
- `method::Mfoil`: The Mfoil method configuration object (unused here but kept for interface consistency).
- `panel_geometry`: Named tuple or mutable structure to be populated with panel data.
- `coordinates::Matrix{TF}`: Nx2 matrix containing x and y coordinates of the airfoil nodes.

# Returns
- The updated `panel_geometry` structure filled with panel edges, vectors, lengths, and node coordinates.
"""
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
