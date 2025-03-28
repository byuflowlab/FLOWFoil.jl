using FLOWMath, Plots

include("geometry_manipulations/utils.jl")

"""
    NACA_65(
        series_number,
        cli;
        N = 161,
        slope_tolerance = 0.000001,
        a = 1.0,
        smoothing = false
    )

This computes the x and y coordinates for any NACA 65 series which can be referenced with a thickness form. It ouptuts the x and y coordinates starting from the trailing edge moving clockwise

# Arguments:
- `series_number::String` : The digits of the airfoil series family, after the 65 (including secondary dashes). For example the NACA 65-3-018 series number would be input as "3-018". Important note, for one reason or another, all airfoils in this family tend to have a 0 after the 2nd dash. For example the NACA 65-3-818 airfoil using the NACA 65-3-018 series when computing airfoil thickness. For this input, if one wanted the NACA 65-3-818 airfoil, they would still put "3-018" for this input. For more information see the NACA Summary of Airfoil Data (Abbott et al, 1945).
- `cli::TF` : Design lift coefficient in tenths of chord. Usually first number after the 2nd dash (ie NACA 65-3-818 would input 0.8 for cli)

# Keyword Arguments:
- `N::Int = 161` : Number of airfoil points (not counting leading edge). Total number of points would be N + num_leading_edge_points
- `slope_tolerance::Float = 0.000001` : The closer this is to 0, the more accurate the derivatives (used in computing airfoil points) should be.
- `a::TF` : Type of mean line used ine chord (if it doesn't say next to the airfoil, it is implied that a = 1.0) - (Abbott et al, 1945)
- `smoothing::Boolean = false` : Creates some extra smoothing on the leading edge - it really only helps for a high number of points (ie ~> 250)

# Returns:
- `x_ouput::Vector{Float}` : Vector of x values for the airfoil that start at the trailing edge and move clockwise.
- `y_output::Vector{Float}` : Vector of y values for the airfoil that start at the trailing edge and move clockwise.
"""
function NACA_65(
    series_number, #series number is a exact string of all the digits after the 65 (including dashes). For example NACA 65-3-018 (Summary of Airfoil Data 1945) would be input as "3-018" 
    cli;
    N = 161,
    slope_tolerance = 0.000001,
    a = 1.0,
    smoothing = false
    )
    # get x coordinates
    N = Int(ceil(N / 2))
    x = [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]

    TF = eltype(x)

    #grabs leading edge radius from the thickness form tables and changes the percentage to a value for computing as well as the center of the circle
    nil, leading_edge_radius = thickness_coordinates(series_number, 0.5)
    leading_edge_radius = leading_edge_radius*(10^-2)
    cos_theta_initial = (0.00501 - 0.005) / sqrt((camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2) #see Summary of Airfoil Data (1945) page 4 for how they compute the leading edge
    leading_edge_circle_center_x = leading_edge_radius*cos_theta_initial#see Summary of Airfoil Data (1945) page 4 for how they compute the leading edge
    #leading_edge_circle_center_y = leading_edge_radius*(camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005)) / sqrt((camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2)

    #transform the cosine spacing to account for the leading edge circle
    x = linear_transform((0, 1), (-leading_edge_radius*cos_theta_initial + leading_edge_circle_center_x, 1), x)

    #define slope variables (theta is the mean line (camber line) slope)
    tan_theta = Vector{TF}(undef, N) .= 0.0
    sin_theta = similar(tan_theta) .= 0.0
    cos_theta = similar(tan_theta) .= 0.0

    #define camber and thickness vectors
    y_c = similar(tan_theta) .= 0.0
    y_t = similar(tan_theta) .= 0.0

    #define output vectors
    y_upper = similar(x) .= 0.0
    y_lower = similar(x) .= 0.0
    x_upper = similar(x) .= 0.0
    x_lower = similar(x) .= 0.0

    
    #compute slope variables
    for i = 2:N-1
        dx_plus = x[i] + slope_tolerance
        dx_minus = x[i] - slope_tolerance
        if x[i] < 0.005 #see Summary of Airfoil Data (1945) page 4
            tan_theta[i] = (camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005)) / (0.00501 - 0.005)
            sin_theta[i] = (camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005)) / sqrt((camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2)
            cos_theta[i] = (0.00501 - 0.005) / sqrt((camber_line(a, cli, 0.00501) - camber_line(a,cli, 0.005))^2 + (0.00501 - 0.005)^2)
        else
            tan_theta[i] = (camber_line(a, cli, dx_plus) - camber_line(a,cli, dx_minus)) / (dx_plus - dx_minus)
            sin_theta[i] = (camber_line(a, cli, dx_plus) - camber_line(a,cli, dx_minus)) / sqrt((camber_line(a, cli, dx_plus) - camber_line(a,cli, dx_minus))^2 + (dx_plus - dx_minus)^2)
            cos_theta[i] = (dx_plus - dx_minus) / sqrt((camber_line(a, cli, dx_plus) - camber_line(a,cli, dx_minus))^2 + (dx_plus - dx_minus)^2)
        end
    end

    if series_number == "-010" #special case exception which uses tabulated data for the slope and camber line function
        ord_1 = [0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25, 30, 35, 40, 45,50, 55, 60, 65, 70, 75, 80, 85, 90, 94, 100].*10^(-2)
        ord_2 = [0.0, 0.250, 0.350, 0.535, 0.930, 1.580, 2.120, 2.585, 3.365, 3.980, 4.475, 4.860, 5.150, 5.355, 5.475, 5.515, 5.475, 5.355, 5.150, 4.860, 4.475, 3.980, 3.365, 2.585, 1.580, 0.0].*10^(-2)
        ord_3 = [0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25, 30, 35, 40, 45,50, 55, 60, 65, 70, 75, 80, 85, 90, 94].*10^(-2)
        slope = [0.42120, 0.38875, 0.34770, 0.29155, 0.23430, 0.19995, 0.17485, 0.13805, 0.11030, 0.08745, 0.06745, 0.04925, 0.03225, 0.01595, 0.0, -0.01595, -0.03225, -0.04925, -0.06745, -0.08745, -0.11030, -0.13805, -0.17485, -0.23430]
        for i = 2:N-1
            if x[i] < 0.005
                y_c[i] = FLOWMath.akima(ord_1, ord_2, 0.005)
                theta = FLOWMath.akima(ord_3, slope, 0.005)
                tan_theta[i] = tan(theta)
                sin_theta[i] = sin(theta)
                cos_theta[i] = cos(theta)
            else
                y_c[i] = FLOWMath.akima(ord_1, ord_2, x[i])
                theta = FLOWMath.akima(ord_3, slope, x[i])
                tan_theta[i] = tan(theta)
                sin_theta[i] = sin(theta)
                cos_theta[i] = cos(theta)
            end
        end
    end

    #compute camber_line for regular cases
    if series_number != "-010"
        for i = 2:N-1
            if x[i] > 0.0
                y_c[i] = camber_line(a, cli, x[i])
            end
        end
    end

    #compute the thickness distributions
    for i = 2:N-1
        if x[i] > 0.0
            y_t[i], nil = thickness_coordinates(series_number, x[i])
        end
    end
    transition_index = 1
    for i = 1:N
        if x[i] < 0.005
            a = sqrt(leading_edge_radius^2 - ((leading_edge_circle_center_x - x[i]) / cos_theta[2])^2)
            x_upper[i] = x[i] - a*sin_theta[i]
            x_lower[i] = x[i] + a*sin_theta[i]
            y_upper_1 = x[i]*tan_theta[i] + a*cos_theta[i]
            y_upper_2 = y_c[i] + y_t[i]*cos_theta[i]
            y_lower_1 = x[i]*tan_theta[i] - a*cos_theta[i]
            y_lower_2 = y_c[i] - y_t[i]*cos_theta[i]
            if y_upper_1 >= y_upper_2
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
            x_upper[i] = x[i] - y_t[i]*sin_theta[i]
            y_upper[i] = y_c[i] + y_t[i]*cos_theta[i]
            x_lower[i] = x[i] + y_t[i]*sin_theta[i]
            y_lower[i] = y_c[i] - y_t[i]*cos_theta[i]
        end
    end
    #gets rid of repeat point
    popfirst!(x_lower)
    popfirst!(y_lower)

    #perform smoothing
    if smoothing == true
        smoothing_interval_top = 3 #number of points that will be linearly interpolated as the leading edge circle transitions to the rest of the airfoil - this should ideally be as small as possible
        i1 = transition_index - round(Int, smoothing_interval_top / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_upper = (y_upper[i2] - y_upper[i1]) / (x_upper[i2] - x_upper[i1])
        for i = i1:i2 
            y_upper[i] = slope_upper*(x_upper[i] - x_upper[i2]) + y_upper[i2]
        end
        smoothing_interval_bottom = 3
        i1 = transition_index - round(Int, smoothing_interval_bottom / 2)
        i2 = transition_index + (smoothing_interval_top - round(Int, smoothing_interval_top / 2))
        slope_lower = (y_lower[i2] - y_lower[i1]) / (x_lower[i2] - x_lower[i1])
        for i = i1:i2
            y_lower[i] = slope_lower*(x_lower[i] - x_lower[i2]) + y_lower[i2]
        end
    end
    #combine vectors into single output starting from the leading edge and moving clockwise
    x_output = [reverse(x_lower); x_upper]
    y_output = [reverse(y_lower); y_upper]
    return x_output, y_output
end

function camber_line(
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

function thickness_coordinates(
    series_number,
    xpt #value of x that the desired y_t is
)
    y_t = 0.0
    leading_edge_radius = 0.0
    if series_number == "3-018"
        x = [0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25,30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100].*10^(-2)
        y = [0.0, 1.324, 1.599, 2.004, 2.728, 3.831, 4.701, 5.424, 6.568, 7.434, 8.093, 8.568, 8.868, 8.990, 8.916, 8.593, 8.045, 7.317, 6.450, 5.486, 4.456, 3.390, 2.325, 1.324, 0.492, 0.0].*10^(-2) 
        leading_edge_radius = 1.92
        y_t = FLOWMath.akima(x,y, xpt)
    elseif  series_number == "-010"
        x = [0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100].*10^(-2)
        y = [0.0, 0.752, 0.890, 1.124, 1.571, 2.222, 2.709, 3.111, 3.746, 4.218, 4.570, 4.824, 4.982, 5.057, 5.029, 4.870, 4.570, 4.151, 3.627, 3.038, 2.451, 1.847, 1.251, 0.749, 0.354, 0.150].*10^(-2)
        leading_edge_radius = 0.666
        y_t = FLOWMath.akima(x,y, xpt)
    end
    return y_t, leading_edge_radius
end