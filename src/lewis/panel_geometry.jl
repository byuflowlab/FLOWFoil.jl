"""
    generate_panel_geometry(method::Lewis, coordinates::AbstractVector)

Generates panel geometries for multiple airfoil surfaces using the Lewis method.

# Arguments
- `method::Lewis`: Configuration or marker for the Lewis panel method.
- `coordinates::AbstractVector{<:Matrix{<:Real}}`: A vector of 2D coordinate matrices, each of size (N+1, 2), defining the panel endpoints for each airfoil surface.

# Returns
- `Vector{NamedTuple}`: A vector of named tuples, each representing the panel geometry for one airfoil. Each tuple includes:
  - `npanels::Int`: Number of panels.
  - `panel_center::Matrix`: Midpoints of each panel.
  - `panel_length::Vector`: Lengths of each panel.
  - `panel_normal::Matrix`: Unit normal vectors to each panel.
  - `panel_curvature::Vector`: Curvature of each panel.
  - `panel_angle::Vector`: Angle of each panel with respect to the x-axis (in radians).
"""
function generate_panel_geometry(method::Lewis, coordinates)
    #broadcast for multiple airfoils
    return [generate_panel_geometry(method, c) for c in coordinates]
end

"""
    generate_panel_geometry(method::Lewis, coordinates::Matrix{<:Real})

Generates panel geometry for a single airfoil surface using the Lewis method.

# Arguments
- `method::Lewis`: Configuration or marker for the Lewis panel method.
- `coordinates::Matrix{<:Real}`: A matrix of size (N+1, 2), where each row is an (x, r) coordinate of a panel corner. The second column must be non-negative (r ≥ 0), corresponding to the radial direction in axisymmetric flow.

# Returns
- `NamedTuple`: A named tuple containing panel geometry fields:
  - `npanels::Int`: Number of panels.
  - `panel_center::Matrix`: Midpoints of each panel.
  - `panel_length::Vector`: Lengths of each panel.
  - `panel_normal::Matrix`: Unit normal vectors for each panel.
  - `panel_curvature::Vector`: Curvature at each panel.
  - `panel_angle::Vector`: Panel orientation angles (in radians).
"""
function generate_panel_geometry(method::Lewis, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), coordinates[:, 2])

    npanels = size(coordinates, 1) - 1

    panel_geometry = (;
        npanels=npanels,
        # - Initialize Outputs - #
        panel_center=zeros(TF, npanels, 2),
        panel_length=zeros(TF, npanels),
        panel_normal=zeros(TF, npanels, 2),
        panel_curvature=zeros(TF, npanels),
        panel_angle=zeros(TF, npanels),
    )

    return generate_panel_geometry!(method, panel_geometry, coordinates)
end

"""
    generate_panel_geometry!(::Lewis, panel_geometry, coordinates::Matrix{<:Real})

In-place computation of panel geometry fields for a single airfoil surface using the Lewis method. This function fills in the provided `panel_geometry` named tuple based on panel corner coordinates.

# Arguments
- `method::Lewis`: Configuration or marker for the Lewis panel method.
- `panel_geometry::NamedTuple`: Preallocated named tuple with fields:
  - `npanels::Int`
  - `panel_center::Matrix`
  - `panel_length::Vector`
  - `panel_normal::Matrix`
  - `panel_curvature::Vector`
  - `panel_angle::Vector`
- `coordinates::Matrix{<:Real}`: (N+1, 2) matrix of panel endpoint coordinates (x, r), where r must be non-negative.

# Returns
- `NamedTuple`: The updated `panel_geometry` tuple with computed values.
"""
function generate_panel_geometry!(
    ::Lewis, panel_geometry, coordinates::Matrix{TF}
) where {TF}

    ### --- SETUP --- ###

    (; npanels, panel_center, panel_length, panel_normal, panel_curvature, panel_angle) =
        panel_geometry

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), r)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_geometry.panel_center[i, :] .= [
            0.5 * (x[i] + x[i + 1])
            0.5 * (r[i] + r[i + 1])
        ]

        # Calculate panel_geometry.panel length
        panel_vector, panel_geometry.panel_length[i] = get_d([x[i] r[i]; x[i + 1] r[i + 1]])

        # Calculate panel_geometry.panel unit normal
        panel_geometry.panel_normal[i, :] = get_panel_normal(
            panel_vector, panel_geometry.panel_length[i]
        )

        # - Calculate panel_geometry.panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        if panel_vector[1] == 0.0
            beta = pi / 2.0
        else
            beta = atan(panel_vector[2] / panel_vector[1])
        end

        # Apply corrections as needed based on orientation of panel_geometry.panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panel_geometry.panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_geometry.panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panel_geometry.panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_geometry.panel_angle[i] = beta + pi
        else
            panel_geometry.panel_angle[i] = beta
        end
    end

    # - Calculate Panel Curvature - #
    # - If not the end panel_geometry, calculate non-zero curvatures - #
    # NOTE: end panel_geometry are trailing edges, which are assumed to have zero curvature.
    for i in 2:(npanels - 1)
        panel_geometry.panel_curvature[i] =
            (panel_geometry.panel_angle[i + 1] - panel_geometry.panel_angle[i - 1]) / 8.0 /
            pi
    end

    return panel_geometry
end
