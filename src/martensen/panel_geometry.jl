"""
    generate_panel_geometry(method::Martensen, coordinates)

Broadcasts `generate_panel_geometry!` over multiple sets of coordinates to generate panel geometries for multiple airfoils.

# Arguments
- `method::Martensen`: The panel method type.
- `coordinates`: A collection (e.g., vector) of coordinate matrices, where each matrix defines the vertices of one airfoil.

# Returns
- A collection of updated `panel_geometry` objects corresponding to each set of input coordinates.
"""
function generate_panel_geometry(method::Martensen, coordinates)
    #broadcast for multiple airfoils
    return generate_panel_geometry.(Ref(method), coordinates)
end

function generate_panel_geometry(method::Martensen, coordinates::AbstractMatrix{TF}) where {TF}

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    npanels = size(coordinates, 1) - 1

    # - Initialize Outputs - #

    panel_geometry = (;
        npanels=npanels,
        # - Initialize Outputs - #
        panel_center=zeros(TF, npanels, 2),
        panel_length=zeros(TF, npanels),
        panel_normal=zeros(TF, npanels, 2),
        panel_angle=zeros(TF, npanels),
        delta_angle=zeros(TF, npanels),
        sine_vector=zeros(TF, npanels),
        cosine_vector=zeros(TF, npanels),
    )

    return generate_panel_geometry!(method, panel_geometry, coordinates)
end

"""
    generate_panel_geometry!(method::Martensen, panel_geometry, coordinates::Matrix{TF}) where {TF}

Computes and updates the geometric properties of panels defined by given coordinates for the Martensen method.

# Arguments
- `method::Martensen`: The panel method type (used here for dispatch, but not directly referenced).
- `panel_geometry`: A named tuple or struct containing fields for panel properties, which will be updated in-place. Expected fields:
    - `npanels`: Number of panels
    - `panel_center`: Matrix to hold panel center coordinates (control points)
    - `panel_length`: Vector to hold lengths of each panel
    - `panel_normal`: Matrix to hold unit normal vectors for each panel
    - `delta_angle`: Vector to hold the change in panel angle (between neighboring panels)
    - `panel_angle`: Vector to hold the angle of each panel relative to the x-axis
    - `sine_vector`: Vector of sines of panel orientation angles
    - `cosine_vector`: Vector of cosines of panel orientation angles
- `coordinates::Matrix{TF}`: An (N+1)×2 matrix of x and y coordinates defining the panel vertices. Each consecutive pair defines a panel.

# Returns
- The updated `panel_geometry` object with all geometric properties filled.
"""
function generate_panel_geometry!(
    method::Martensen, panel_geometry, coordinates::Matrix{TF}
) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Unpack Tuple
    (;
        npanels,
        panel_center,
        panel_length,
        panel_normal,
        delta_angle,
        panel_angle,
        sine_vector,
        cosine_vector,
    ) = panel_geometry

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]])
        # panel_length[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)
        sine_vector[i] = (y[i + 1] - y[i]) / panel_length[i]
        cosine_vector[i] = (x[i + 1] - x[i]) / panel_length[i]

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # compute the panel angle using Lewis' method of doing so (See program 1.1 data preparation in Lewis 1991)
        abscosine = abs(cosine_vector[i])
        tol = sqrt(eps()) # lewis uses 0.000001 #term used to compute the panel angle (slope)
        if abscosine > tol
            t = atan(sine_vector[i] / cosine_vector[i])
        end
        if abscosine <= tol
            panel_angle[i] = (sine_vector[i] / abs(sine_vector[i])) * pi / 2
        end
        if cosine_vector[i] > tol
            panel_angle[i] = t
        end
        if cosine_vector[i] < -tol
            panel_angle[i] = t - pi
        end
    end

    for i in 1:npanels
        if i == 1 || i == npanels
            delta_angle[i] = (panel_angle[i]) / 2.0
        else
            delta_angle[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 2.0
        end
    end
    # - Return Panel Object - #
    return panel_geometry
end
