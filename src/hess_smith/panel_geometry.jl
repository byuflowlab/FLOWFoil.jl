#=
panel_geometry = (; control_points, panel_lengths, panel_angles, ??)
copy and paste from martensen and change.
=#
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
        sine_vector=zeros(TF, npanels),
        cosine_vector=zeros(TF, npanels),
    )

    return generate_panel_geometry!(method, panel_geometry, coordinates)
end

function generate_panel_geometry!(
    method::HessSmith, panel_geometry, coordinates::Matrix{TF}
) where {TF}

    # Unpack Tuple
    (;
        npanels,
        x,
        y,
        panel_center,
        panel_length,
        panel_angle,
        sine_vector,
        cosine_vector,
    ) = panel_geometry

    x = coordinates[:, 1]
    y = coordinates[:, 2]

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]])
        sine_vector[i] = (y[i + 1] - y[i]) / panel_length[i]
        cosine_vector[i] = (x[i + 1] - x[i]) / panel_length[i]
        panel_angle[i] = asind(sine_vector[i])
    end

    # - Return Panel Object - #
    return panel_geometry
end
