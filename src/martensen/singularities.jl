"""
    calculate_periodic_vortex_influence(paneli, panelj)

Cacluate the influence of a periodic vortex at panel j onto panel i.

# Arguments:
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

# Returns:
- `aij::Float` : Influence of vortex strength at panel j onto panel i.
"""
function calculate_periodic_vortex_influence(
    paneli, panelj, system_geometry, i, j, cascade_parameters
)

    # - Rename for convenience - #
    # control points
    r_x = system_geometry.r_x[i, j]
    r_y = system_geometry.r_y[i, j]

    # cascade stuff
    stagger = cascade_parameters.stagger
    pitch = cascade_parameters.pitch

    # adjusted angles
    sine_panel_angle_i = sin(stagger + paneli.panel_angle[i])
    cosine_panel_angle_i = cos(stagger + paneli.panel_angle[i])

    # - Coefficient - #
    if pitch < eps()
        a = b = k = 0.0
    else
        a = -(r_x * cos(stagger) - r_y * sin(stagger)) * 2.0 * pi / pitch
        b = -(r_x * sin(stagger) + r_y * cos(stagger)) * 2.0 * pi / pitch
        cosha = 0.5 * (exp(a) + 1.0 / exp(a))
        k = 0.5 / pitch / (cosha - cos(b))
    end
    sinha = 0.5 * (exp(a) - 1.0 / exp(a))

    return (sinha * sine_panel_angle_i - sin(b) * cosine_panel_angle_i) *
           k *
           panelj.panel_length[j]
end

function calculate_planar_vortex_influence(
    paneli, panelj, system_geometry, i, j, cascade_parameters
)

    # rename for convenience
    x_i = paneli.panel_center[i]
    x_j = panelj.panel_center[j]
    y_i = paneli.panel_center[i]
    y_j = panelj.panel_center[j]

    # cascade stuff
    stagger = cascade_parameters.stagger

    # adjusted angles
    sine_panel_angle_i = sin(stagger + paneli.panel_angle[i])
    cosine_panel_angle_i = cos(stagger + paneli.panel_angle[i])

    # u = (y_j - y_i) / (r * 2 * pi)
    u = system_geometry.r_y[i, j] / (2.0 * pi * system_geometry.r_squared[i, j])

    # v = -(x_j - x_i) / (r * 2 * pi)
    v = -system_geometry.r_x[i, j] / (2.0 * pi * system_geometry.r_squared[i, j])

    return (u * cosine_panel_angle_i + v * sine_panel_angle_i) * panelj.panel_length[j]
end

function calculate_periodic_self_vortex_influence(panel, i)
    #compute self-inducing coupling coefficients
    # return -0.5 - (2.0 * panel.delta_angle[i] - 2.0 * pi) / (8 * pi)
    return 0.5 - (panel.delta_angle[i] - pi) / (4.0 * pi)
    # return 0.5
end
