"""
    NACA4

# Fields
- `maximum_camber::Float=2.0` : maximum camber in % chord
- `maximum_camber_position::Float=4.0` : x-position of maximum camber point in 1/10 chord
- `maximum_thickness::Float=12.0` : maximum thickness in % chord
- `blunt_trailing_edge::Bool=false` : flag for whether to use the blunt trailing edge NACA 4-series definition or not.
"""
@kwdef struct NACA4{Tc,Tp,Tt,Tb} <: AirfoilGeometry
    maximum_camber::Tc = 2.0
    maximum_camber_position::Tp = 4.0
    maximum_thickness::Tt = 12.0
    blunt_trailing_edge::Tb = false
end

"""
    naca4_thickness(x, maximum_thickness; blunt_trailing_edge=false)

Compute thickness at a given chord-normalized x-position by NACA 4-series thickness equations.

# Arguments
- `x::Float` : x position along chordline, markersize=3, markershape=:squaree
- `maximum_thickness::Float` : Maximum thickness value

# Keyword Arguments
- `blunt_trailing_edge::Bool=false` : Flag whether trailing edge is blunt or not
"""
function naca4_thickness(x, maximum_thickness; blunt_trailing_edge=false)

    # change c2 coefficient base on whether trailing edge is to be blunt or not.
    c2 = blunt_trailing_edge ? 0.3516 : 0.3537

    # return naca4_thickness value
    return 10.0 *
           maximum_thickness *
           (0.2969 * sqrt(x) - 0.1260 * x - c2 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

"""
    naca4_camber(x, maximum_camber, maximum_camber_position)

Compute camber at a given chord-normalized x-position by NACA 4-series camber equations.

# Arguments
- `x::Float` : x position along chordline
- `maximum_camber::Float64` : Maximum camber value
- `maximum_camber_position::Float64` : Position of maximum camber
"""
function naca4_camber(x, maximum_camber, maximum_camber_position)
    if real(maximum_camber) != 0.0 && real(maximum_camber_position) != 0.0
        if x <= real(maximum_camber_position)
            zbar =
                maximum_camber * (2 * maximum_camber_position * x - x^2) /
                maximum_camber_position^2
        else
            zbar =
                maximum_camber *
                (1 - 2 * maximum_camber_position + 2 * maximum_camber_position * x - x^2) /
                (1 - maximum_camber_position)^2
        end
    else
        zbar = 0.0
    end
    return zbar
end

"""
    naca4(parameters::NACA4; N=161, x=nothing, split=false)

Compute x, y airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

# Arguments
- `parameters::NACA4` : NACA 4-series parameters

# Keyword Arguments
- `N::Int=161` : Total number of coordinates to use.  This values should be odd, but if not, the number of points returned will be N-1.
- `x::AbstractArray{Float}` : x coordinates (cosine spaced coordinates used by default)
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function naca4(parameters::NACA4; N=161, x=nothing, split=false)
    return naca4(
        parameters.maximum_camber,
        parameters.maximum_camber_position,
        parameters.maximum_thickness;
        N=N,
        x=x,
        blunt_trailing_edge=parameters.blunt_trailing_edge,
        split=split,
    )
end

"""
    naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_trailing_edge=false, split=false)

Compute x, y airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

# Arguments
- `c::Float` : Maximum camber value (percent of chord)
- `p::Float` : Position along chord (in 10ths of chord) where maximum naca4_camber lies
- `t::Float` : Maximum thickness of airfoil in percent chord

# Keyword Arguments
- `N::Int` : Total number of coordinates to use.  This values should be odd, but if not, the number of points returned will be N-1.
- `x::AbstractArray{Float}` : x-coordinates (cosine spaced coordinates used by default)
- `blunt_trailing_edge::Bool` : Flag whether trailing edge is blunt or not
- `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function naca4(
    c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_trailing_edge=false, split=false
)

    # get x coordinates
    if isnothing(x)
        N = ceil(Int, N / 2)
        x = split_cosine_spacing(N)
    else
        N = ceil(Int, length(x) / 2)
    end

    #naca digits
    maximum_camber = c / 100.0
    maximum_camber_position = p / 10.0
    maximum_thickness = t / 100.0

    #initialize arrays
    TF = promote_type(typeof(c), eltype(x))
    zu = zeros(TF, N) #upper y values
    zl = zeros(TF, N) #lower y values

    #--Calculate y-values--#
    #naca4_thickness distribution
    T = naca4_thickness.(x, Ref(maximum_thickness); blunt_trailing_edge=blunt_trailing_edge)

    #naca4_camber distribution
    zbar = naca4_camber.(x, Ref(maximum_camber), Ref(maximum_camber_position))

    #y-positions at chordwise stations
    zl = @. zbar - T / 2
    zu = @. zbar + T / 2

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_naca4(x,y)

Calculate NACA 4-series parameters based on input x,y coordinates.

# Arguments
- `x::AbstractArray{Float}` : vector of x coordinates for airfoil
- `y::AbstractArray{Float}` : vector of y coordinates for airfoil

# Keyword Arguments
- `blunt_trailing_edge::Bool=false` : Flag whether trailing edge is blunt or not

# Returns
- `parameters::NACA4` : a parameter object of type NACA4.
"""
function determine_naca4(x, y; blunt_trailing_edge=false)

    # model to fit
    function model(x, param)
        _, y = naca4(
            param[1],
            param[2],
            param[3];
            x=x[floor(Int, length(x) / 2)+1:end],
            blunt_trailing_edge=blunt_trailing_edge,
        )
        return y
    end

    # solve for coefficients
    ufit = LsqFit.curve_fit(model, x, y, [1.0, 1.0, 10.0])

    # unpack
    return NACA4(;
        maximum_camber=ufit.param[1],
        maximum_camber_position=ufit.param[2],
        maximum_thickness=ufit.param[3],
    )
end
