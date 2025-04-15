######################################################################
#                                                                    #
#               NACA 65-series Compressor Airfoil                    #
#                                                                    #
######################################################################
"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 1
"""
function scaled_thickness(x; method="scaled")
    if method == "derived"
        tr =
            [
                0.0
                0.772
                0.932
                1.169
                1.574
                2.177
                2.647
                3.040
                3.666
                4.143
                4.503
                4.760
                4.924
                4.996
                4.963
                4.812
                4.530
                4.146
                3.682
                3.156
                2.584
                1.987
                1.385
                0.810
                0.306
                0.0
            ] * 1e-2
        ler = 0.687
    elseif method == "scaled"
        tr =
            [
                0.0
                0.752
                0.890
                1.124
                1.571
                2.222
                2.709
                3.111
                3.746
                4.218
                4.570
                4.824
                4.982
                5.057
                5.029
                4.870
                4.570
                4.151
                3.627
                3.038
                2.451
                1.847
                1.251
                0.749
                0.354
                0.150
            ] * 1e-2
        ler = 0.666
    else
        @error "no method $method, please choose scaled or derived."
    end
    tx =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2

    return FLOWMath.akima(tx, tr, x)
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function scaled_camber(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    a1y =
        [
            0.0
            0.25
            0.35
            0.535
            0.93
            1.580
            2.120
            2.585
            3.365
            3.98
            4.475
            4.86
            5.15
            5.355
            5.475
            5.515
            5.475
            5.355
            5.15
            4.86
            4.475
            3.98
            3.365
            2.585
            1.58
            0.0
        ] * 1e-2

    a1fine = FLOWMath.akima(a1x, a1y, x)

    return a1fine * clo
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function scaled_slope(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    dydx = [
        0.4212
        0.38875
        0.3477
        0.29155
        0.2343
        0.19995
        0.17485
        0.13805
        0.11030
        0.08745
        0.06745
        0.04925
        0.03225
        0.01595
        0.0
        -0.01595
        -0.03225
        -0.04925
        -0.06745
        -0.08745
        -0.11030
        -0.13805
        -0.17485
        -0.23430
    ]

    dydxfine = FLOWMath.akima(a1x[2:(end - 1)], dydx, x)

    return dydxfine * clo
end

"""
    computed_camber_line(
    a,
    cli,
    x
)

This computes the corresponding camber for a given x value for the NACA 65 series. 

# Arguments:
- `a::TF` : Mean-line designation, fraction of chord from leading edge over which design load is uniform. Should be listed in the airfoil descirption
- `cli::TF` : Design lift coefficient in tenths of chord. Usually first number after the 2nd dash (ie NACA 65-3-818 would input 0.8 for cli)
- `x::TF` : x value (note that x must be normalized by the chord length)

# Returns:
- `yc::Float` : Amount of camber y_c for a given x value on the airfoil
"""
function computed_camber_line(
    a,
    cli,
    x
)
    if a == 1.0 #There is a discountinuity when a = 1 that must be removed.
        a = 0.9999
    end
    g = (- 1 / (1 - a))*(
        (a^2)*(0.5*log(a) - 0.25) + 0.25
    )
    h = g + (1 / (1 - a))*(
        0.5*((1-a)^2)*log(1 - a) - 0.25*(1 - a)^2
    )
    y_c = (cli / (2*pi*(a + 1)))*(
        (1 / (1 - a))*(0.5*((a-x)^2)*log(abs(a-x)) - 0.5*((1 - x)^2)*log(1 - x) + 0.25*(1-x)^2 - 0.25*(a - x)^2)
        - x*log(x) + g - h*x
    )
    return y_c
end

"""
    slopes65(
    x,
    clo,
    a;
    slope_tolerance = 0.000001
)

This computes the corresponding slopes of the mean line for a given x value on the airfoil

# Arguments:
- `x::TF` : x value (note that x must be normalized by the chord length)
- `clo::TF` : Design lift coefficient in tenths of chord. Usually first number after the 2nd dash (ie NACA 65-3-818 would input 0.8 for cli)
- `a::TF` : Mean-line designation, fraction of chord from leading edge over which design load is uniform. Should be listed in the airfoil descirption

# Keyword Arguments:
- `slope_tolerance::Float = 0.000001` : step size taken for the derivative

# Returns:
- `tan_theta::Float` : Tangent of the mean line slope
- `sin_theta::Float` : Sine of the mean line slope
- `cos_theta::Float` : Cosine of the mean line slope
"""
function slopes65(
    x,
    clo,
    a;
    slope_tolerance = 0.000001
)
    cli = clo
    dx_plus = x + slope_tolerance
    dx_minus = x - slope_tolerance
    tan_theta = 0.0
    sin_theta = 0.0
    cos_theta = 0.0
    if x < 0.005 #see Summary of Airfoil Data (1945) page 4
        tan_theta = (computed_camber_line(a, cli, 0.00501) - computed_camber_line(a,cli, 0.005)) / (0.00501 - 0.005)
        sin_theta = (computed_camber_line(a, cli, 0.00501) - computed_camber_line(a,cli, 0.005)) / sqrt((computed_camber_line(a, cli, 0.00501) - computed_camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2)
        cos_theta = (0.00501 - 0.005) / sqrt((computed_camber_line(a, cli, 0.00501) - computed_camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2)
    else
        tan_theta = (computed_camber_line(a, cli, dx_plus) - computed_camber_line(a,cli, dx_minus)) / (dx_plus - dx_minus)
        sin_theta = (computed_camber_line(a, cli, dx_plus) - computed_camber_line(a,cli, dx_minus)) / sqrt((computed_camber_line(a, cli, dx_plus) - computed_camber_line(a,cli, dx_minus))^2 + (dx_plus - dx_minus)^2)
        cos_theta = (dx_plus - dx_minus) / sqrt((computed_camber_line(a, cli, dx_plus) - computed_camber_line(a,cli, dx_minus))^2 + (dx_plus - dx_minus)^2)
    end
    return tan_theta, sin_theta, cos_theta
end

"""
    thickness65(
    series_number,
    xpt #value of x that the desired y_t is
)

This computes the thicknes coordinate and leading edge radius from tabulated thickness forms

# Arguments:
- `series_number::String` : Airfoil series number (see naca65 doc string for more details)
- `xpt::TF` : Desired x value to compute the thickness value

# Returns:
- `y_t::Float` : thickness value at specified x value
- `leading_edge_radius::Float` : Tabulated leading edge thickness - make sure it is in x/c units!
"""
function thickness65(
    series_number,
    xpt #value of x that the desired y_t is
)
    y_t = 0.0
    leading_edge_radius = 0.0
    if series_number == "3-018"
        x = [0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100].*10^(-2)
        y = [0.0, 1.324, 1.599, 2.004, 2.728, 3.831, 4.701, 5.424, 6.568, 7.434, 8.093, 8.568, 8.868, 8.990, 8.916, 8.593, 8.045, 7.317, 6.450, 5.486, 4.456, 3.390, 2.325, 1.324, 0.492, 0.0].*10^(-2) 
        leading_edge_radius = 1.92*(10^(-2))
        y_t = FLOWMath.akima(x,y, xpt)
    end
    return y_t, leading_edge_radius
end

"""
Assumes x in non-dimensional range [0.0,1.0]


Description from NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds:"

The 65-series compressor blade family is formed by combining a basic thickness form with cambered mean lines.
The basic thickness form used is the NACA 65(216)-010 thickness form with the ordinates increased by 0.0015 times the chordwise stations to provide slightly, increased thickness toward the trailing edge.
In the scaled case, it was not derived for 10-percent thickness but was scaled down from the NACA 65,2-016 airfoil.
The scaling procedure gives the best results whep it is restricted to maximum thickness changes of a few percent.
The NACA 65-010 basic thickness has also been derived.
These thickness forms differ slightly but are considered to be interchangeable.

The basic mean line used is the a=1.0 mean line.
The amount of camber is for the design lift coefficient for the isolated airfoil with cl_o of 1.0.
Both ordinates and slopes are scaled directly to obtain other cambers.
Cambered blade sections are obtained by applying the thickness perpendicular to the mean line at stations laid out along the chord line.
In the designation the camber is given by the first number after the dash in tenths of cl_o.
For example, the NACA 65-810 and NACA 65-(12)10 blade sections are cambered for cl_o = 0.8 and cl_o = 1.2, respectively.

# Arguments:
- `clo::TF` : Design lift coefficient in tenths of chord. Usually first number after the 2nd dash (ie NACA 65-3-818 would input 0.8 for clo)
- `a::TF` : Mean-line designation, fraction of chord from leading edge over which design load is uniform.
- `series_number::String` : digits of the airfoil series family with clo of 0 (ie NACA 65-3-818 would be "3-010")

# Keyword Arguments:
- `x::Vector{TF} = nothing` : input x values if specificing x values manually
- `split::Boolean = false` : if true, then the output will be split between top and bottom coordinates
- `extra_blending::Boolean = false` : If desired number of points is large (> 300ish) then set to true and it will add some extra blending if desired. Note: This is generally not needed!

# ouptuts:
- `x::Vector{TF}` : x coordinates
- `y::Vector{TF}` : y coordinates
"""
function naca65(clo, a, series_number; N=161, x=nothing, split=false, extra_blending=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = cosine_spacing(N)
    end
    nil, leading_edge_radius = thickness65(series_number, 0.5)
    tan_theta_initial, sin_theta_initial, cos_theta_initial = slopes65(0.001, clo, a)
    leading_edge_circle_center_x = leading_edge_radius * cos_theta_initial
    x = linear_transform(
        (0, 1),
        (-leading_edge_radius * cos_theta_initial + leading_edge_circle_center_x, 1),
        x,
    )
    y_upper = similar(x) .= 0.0
    y_lower = similar(x) .= 0.0
    x_upper = similar(x) .= 0.0
    x_lower = similar(x) .= 0.0
    cos_theta = similar(x) .= 0.0
    sin_theta = similar(x) .= 0.0
    tan_theta = similar(x) .= 0.0
    c = similar(x) .= 0.0
    t = similar(x) .= 0.0
    transition_index = 1 #this is used for extra blending - it is the index when the leading edge circle transitions to the rest of the airfoil

    for i in 1:N
        if x[i] < 0.005
            tan_theta[i] = tan_theta_initial
            sin_theta[i] = sin_theta_initial
            cos_theta[i] = cos_theta_initial
            c[i] = computed_camber_line(a, clo, x[i])
            t[i], nil = thickness65(series_number, x[i])
            rad_a =
                leading_edge_radius^2 -
                ((leading_edge_circle_center_x - x[i]) / cos_theta_initial)^2
            alpha = sqrt(rad_a)
            x_upper[i] = x[i] - alpha * sin_theta[i]
            x_lower[i] = x[i] + alpha * sin_theta[i]
            y_upper_1 = x[i] * tan_theta[i] + alpha * cos_theta[i]
            y_upper_2 = c[i] + t[i] * cos_theta[i]
            y_lower_1 = x[i] * tan_theta[i] - alpha * cos_theta[i]
            y_lower_2 = c[i] - t[i] * cos_theta[i]
            if y_upper_1 >= y_upper_2
                y_upper[i] = y_upper_1
            elseif x[i] < 0.005 / 4
                y_upper[i] = y_upper_1
            else
                y_upper[i] = y_upper_2
            end
            if y_lower_1 <= y_lower_2
                y_lower[i] = y_lower_1
            else
                y_lower[i] = y_lower_2
            end
            transition_index = i
        else
            if x[i] != 1.0
                tan_theta[i], sin_theta[i], cos_theta[i] = slopes65(x[i], clo, a)
                c[i] = computed_camber_line(a, clo, x[i])
                t[i], nil = thickness65(series_number, x[i])
            else
                c[i] = 0.0
                t[i] = 0.0
            end
            x_upper[i] = x[i] - t[i] * sin_theta[i]
            y_upper[i] = c[i] + t[i] * cos_theta[i]
            x_lower[i] = x[i] + t[i] * sin_theta[i]
            y_lower[i] = c[i] - t[i] * cos_theta[i]
        end
    end
    #perform extra blending at the leading edge - good for many points aka >300ish points
    if extra_blending == true
        smoothing_interval_top = 3 #number of points that will be linearly interpolated as the leading edge circle transitions to the rest of the airfoil - this should ideally be as small as possible
        i1 = transition_index - round(Int, smoothing_interval_top / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_upper = (y_upper[i2] - y_upper[i1]) / (x_upper[i2] - x_upper[i1])
        for i = i1:i2 
            y_upper[i] = slope_upper*(x_upper[i] - x_upper[i2]) + y_upper[i2]
        end
        smoothing_interval_bottom = 10 #feel free to adjust this value to get a better blend - each airfoil you kinda have to play around with when you add lots of points
        i1 = transition_index - round(Int, smoothing_interval_bottom / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_lower = (y_lower[i2] - y_lower[i1]) / (x_lower[i2] - x_lower[i1])
        for i = i1:i2
            y_lower[i] = slope_lower*(x_lower[i] - x_lower[i2]) + y_lower[i2]
        end
    end
    if split
        return reverse(x_lower), x_upper, reverse(y_lower), y_upper
    else
        return [reverse(x_lower[2:end]); x_upper[2:end]], [reverse(y_lower[2:end]); y_upper[2:end]]
    end
end

"""
Assumes x in non-dimensional range [0.0,1.0]


Description from NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds:"

The 65-series compressor blade family is formed by combining a basic thickness form with cambered mean lines.
The basic thickness form used is the NACA 65(216)-010 thickness form with the ordinates increased by 0.0015 times the chordwise stations to provide slightly, increased thickness toward the trailing edge.
In the scaled case, it was not derived for 10-percent thickness but was scaled down from the NACA 65,2-016 airfoil.
The scaling procedure gives the best results whep it is restricted to maximum thickness changes of a few percent.
The NACA 65-010 basic thickness has also been derived.
These thickness forms differ slightly but are considered to be interchangeable.

The basic mean line used is the a=1.0 mean line.
The amount of camber is for the design lift coefficient for the isolated airfoil with cl_o of 1.0.
Both ordinates and slopes are scaled directly to obtain other cambers.
Cambered blade sections are obtained by applying the thickness perpendicular to the mean line at stations laid out along the chord line.
In the designation the camber is given by the first number after the dash in tenths of cl_o.
For example, the NACA 65-810 and NACA 65-(12)10 blade sections are cambered for cl_o = 0.8 and cl_o = 1.2, respectively.
"""
function naca65_scaled(clo; N=161, x=nothing, split=false, extra_blending = false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = cosine_spacing(N)
    end

    leading_edge_radius = 0.666*10^(-2)
    theta_initial = scaled_slope(clo, 0.005)
    cos_theta_initial = cos(theta_initial)
    leading_edge_circle_center_x = leading_edge_radius * cos_theta_initial
    x = linear_transform(
        (0, 1),
        (-leading_edge_radius * cos_theta_initial + leading_edge_circle_center_x, 1),
        x,
    )

    #define output vectors
    y_upper = similar(x) .= 0.0
    y_lower = similar(x) .= 0.0
    x_upper = similar(x) .= 0.0
    x_lower = similar(x) .= 0.0

    cos_theta = similar(x) .= 0.0
    sin_theta = similar(x) .= 0.0
    tan_theta = similar(x) .= 0.0
    c = similar(x) .= 0.0
    t = similar(x) .= 0.0
    transition_index = 1

    for i in 1:N
        if x[i] < 0.005
            c[i] = scaled_camber(clo, 0.005)
            t[i] = scaled_thickness(x[i])
            tan_theta[i] = tan(theta_initial)
            sin_theta[i] = sin(theta_initial)
            cos_theta[i] = cos(theta_initial)
            rad_a =
                leading_edge_radius^2 -
                ((leading_edge_circle_center_x - x[i]) / cos_theta[i])^2
            a = sqrt(rad_a)
            x_upper[i] = x[i] - a * sin_theta[i]
            x_lower[i] = x[i] + a * sin_theta[i]
            y_upper_1 = x[i] * tan_theta[i] + a * cos_theta[i]
            y_upper_2 = c[i] + t[i] * cos_theta[i]
            y_lower_1 = x[i] * tan_theta[i] - a * cos_theta[i]
            y_lower_2 = c[i] - t[i] * cos_theta[i]
            if y_upper_1 >= y_upper_2
                y_upper[i] = y_upper_1
            elseif x[i] < 0.005 / 4
                y_upper[i] = y_upper_1
            else
                y_upper[i] = y_upper_2
            end
            if y_lower_1 <= y_lower_2
                y_lower[i] = y_lower_1
            else
                y_lower[i] = y_lower_2
            end
            transition_index = i
        else
            sin_theta[i] = sin(scaled_slope(clo, x[i]))
            cos_theta[i] = cos(scaled_slope(clo, x[i]))
            tan_theta[i] = tan(scaled_slope(clo, x[i]))
            c[i] = scaled_camber(clo, x[i])
            t[i] = scaled_thickness(x[i])
            x_upper[i] = x[i] - t[i] * sin_theta[i]
            y_upper[i] = c[i] + t[i] * cos_theta[i]
            x_lower[i] = x[i] + t[i] * sin_theta[i]
            y_lower[i] = c[i] - t[i] * cos_theta[i]
        end
    end
        #perform smoothing
    if extra_blending == true
        smoothing_interval_top = 3 #number of points that will be linearly interpolated as the leading edge circle transitions to the rest of the airfoil - this should ideally be as small as possible
        i1 = transition_index - round(Int, smoothing_interval_top / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_upper = (y_upper[i2] - y_upper[i1]) / (x_upper[i2] - x_upper[i1])
        for i = i1:i2 
            y_upper[i] = slope_upper*(x_upper[i] - x_upper[i2]) + y_upper[i2]
        end
        smoothing_interval_bottom = 7
        i1 = transition_index - round(Int, smoothing_interval_bottom / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_lower = (y_lower[i2] - y_lower[i1]) / (x_lower[i2] - x_lower[i1])
        for i = i1:i2
            y_lower[i] = slope_lower*(x_lower[i] - x_lower[i2]) + y_lower[i2]
        end
    end

    if split
        return reverse(x_lower), x_upper, reverse(y_lower), y_upper
    else
        return [reverse(x_lower); x_upper[2:end]], [reverse(y_lower); y_upper[2:end]]
    end
end

