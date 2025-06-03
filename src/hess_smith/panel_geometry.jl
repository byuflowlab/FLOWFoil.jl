"""
    generate_panel_geometry(method::HessSmith, coordinates)

Generates panel geometry for an airfoil surface to be used in a Hess-Smith panel method.

# Arguments
- `method::HessSmith`: A marker or configuration type representing the Hess-Smith method.
- `coordinates`: Either a single `Matrix{<:Real}` of size (N+1, 2), where each row is an (x, y) coordinate of a panel corner, or an array of such matrices (to process multiple airfoils).

# Returns
- `panel_geometry::NamedTuple`: A named tuple containing:
  - `npanels::Int`: Number of panels.
  - `x::Vector`: x-coordinates of panel endpoints.
  - `y::Vector`: y-coordinates of panel endpoints.
  - `panel_center::Matrix`: Midpoints of each panel.
  - `panel_length::Vector`: Lengths of each panel.
  - `panel_angle::Vector`: Angles (in degrees) of each panel relative to the x-axis.
  - `panel_vectors::Matrix`: Tangent vectors for each panel.
  - `sine_vector::Vector`: Sine of panel orientation angles.
  - `cosine_vector::Vector`: Cosine of panel orientation angles.
"""
function generate_panel_geometry(method::HessSmith, coordinates)
    #broadcast for multiple airfoils
    return generate_panel_geometry.(Ref(method), coordinates)
end

function generate_panel_geometry(method::HessSmith, coordinates::AbstractMatrix{TF}) where {TF}

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    npanels = size(coordinates, 1) - 1

    # - Initialize Outputs - #

    panel_geometry = (;
        npanels=npanels,
        # - Initialize Outputs - #
        x=zeros(TF, npanels+1),
        y=zeros(TF, npanels+1),
        panel_center=zeros(TF, npanels, 2),
        panel_length=zeros(TF, npanels),
        panel_angle=zeros(TF, npanels),
        panel_vectors=zeros(TF, npanels, 2),
        sine_vector=zeros(TF, npanels),
        cosine_vector=zeros(TF, npanels),
    )

    return generate_panel_geometry!(method, panel_geometry, coordinates)
end

"""
    generate_panel_geometry!(method::HessSmith, panel_geometry, coordinates::Matrix{<:Real})

In-place function that computes and populates the panel geometry for a given set of airfoil surface coordinates,
used in the Hess-Smith panel method.

# Arguments
- `method::HessSmith`: An instance of the Hess-Smith panel method. Currently acts as a marker.
- `panel_geometry::NamedTuple`: A preallocated named tuple containing arrays for storing panel geometry properties.
- `coordinates::Matrix{<:Real}`: A (N+1)×2 matrix where each row is an (x, y) coordinate of a panel corner. Panels are defined between each pair of adjacent rows.

# Computed Values in `panel_geometry`
- `x::Vector`: x-coordinates of panel endpoints.
- `y::Vector`: y-coordinates of panel endpoints.
- `npanels::Int`: Number of panels (equal to number of coordinate points minus one).
- `panel_center::Matrix`: (x, y) coordinates of panel midpoints (control points), size (N, 2).
- `panel_length::Vector`: Length of each panel.
- `panel_angle::Vector`: Orientation angle of each panel in degrees.
- `panel_vectors::Matrix`: Tangent vector of each panel.
- `sine_vector::Vector`: Sine of panel orientation angle.
- `cosine_vector::Vector`: Cosine of panel orientation angle.

# Returns
- `panel_geometry::NamedTuple`: The updated geometry tuple with all computed values filled in.
"""
function generate_panel_geometry!(
    method::HessSmith, panel_geometry, coordinates::Matrix{TF}
) where {TF}

    # Unpack Tuple
    (;
        x,
        y,
        npanels,
        panel_center,
        panel_length,
        panel_angle,
        sine_vector,
        cosine_vector,
    ) = panel_geometry

    for i in 1:npanels+1
        # Define the coordinates
        x[i]=coordinates[i, 1]
        y[i]=coordinates[i, 2]
    end

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        _, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]])
        sine_vector[i] = (y[i + 1] - y[i]) / panel_length[i]
        cosine_vector[i] = (x[i + 1] - x[i]) / panel_length[i]
        panel_angle[i] = asind(sine_vector[i])
        
        #Calculate Panel Vectors
        panel_geometry.panel_vectors[i, :] = [x[i + 1] - x[i]; y[i + 1] - y[i]]
       
    end

    # - Return Panel Object - #
    return panel_geometry
end
