function generate_panels(method::Martensen, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(method), coordinates)
end

function generate_panels(::Martensen, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # - Rename for Convenience - #
    npanels = length(x) - 1

    # - Initialize Outputs - #
    panel_center = zeros(TF, npanels, 2)
    panel_length = zeros(TF, npanels)
    panel_normal = zeros(TF, npanels, 2)
    panel_curvature = zeros(TF, npanels)
    panel_angle = zeros(TF, npanels)
    delta_angle = zeros(TF, npanels)
    sine_vector = zeros(TF. npanels)
    cosine_vector = zeros(TF, npanels)
    ex = 0.000001 #term used to compute the panel angle (slope)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        panel_length[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)
        sine_vector[i] = (y[i + 1] - y[i]) / panel_length[i]
        cosine_vector[i] = (x[i + 1] - x[i]) / panel_length[i]
        
        #panel_vector, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]]) #not sure what this is supposed to be so I commented it out

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        beta = atan(panel_vector[2] / panel_vector[1])

        # Apply corrections as needed based on orientation of panel in coordinate frame.
        #compute the panel angle using Lewis' method of doing so (See program 1.1 data preperation in Lewis 1991)
        abscosine = abs(cosine_vector[i])
        if abscosine > ex
            t = atan(sine_vector[i] / cosine_vector[i])
        end
        if abscosine <= ex
            panel_angle[i] = (sine_vector[i] / abs(sine_vector[i]))*pi / 2
        end
        if cosine_vector[i] > ex
            panel_angle[i] = t
        end
        if cosine_vector[i] < -ex
            panel_angle[i] = t - pi
        end
        #=
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta + pi
        else
            panel_angle[i] = beta
        end
        =#
    end

    for i in 1:npanels
        if i == 1 || i == npanels
            delta_angle[i] = (panel_angle[i]) / 2.0
        else
            delta_angle[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 2.0
        end
    end

    # - Return Panel Object - #
    return ConstantFlatPanel(
        npanels, panel_center, panel_length, panel_normal, panel_angle, delta_angle, sine_vector, cosine_vector
    )
end
