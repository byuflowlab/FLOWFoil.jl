#=
NACA 4-series Parametrization Functions
=#

"""
    thickness(x, maxthick; blunt_TE=false)

Compute thickness at a given chord-normalized x-position by NACA 4-series thickness equations.

**Arguments:**
 - `x::Float` : x position along chordline
 - `maxthick::Float` : Maximum thickness value

**Keyword Arguments:**
 - `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
"""
function thickness(x, maxthick; blunt_TE=false)

    # change c2 coefficient base on whether trailing edge is to be blunt or not.
    c2 = blunt_TE ? 0.3516 : 0.3537

    # return thickness value
    return 10.0 *
           maxthick *
           (0.2969 * sqrt(x) - 0.1260 * x - c2 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

"""
    camber(x, maxcamber, camberpose)

Compute camber at a given chord-normalized x-position by NACA 4-series camber equations.

**Arguments:**
 - `x::Float` : x position along chordline
 - `maxcamber::Float64` : Maximum camber value
 - `camberpose::Float64` : Position of maximum camber
"""
function camber(x, maxcamber, camberpose)
    if real(maxcamber) != 0.0 && real(camberpose) != 0.0
        if x <= real(camberpose)
            ybar = maxcamber * (2 * camberpose * x - x^2) / camberpose^2
        else
            ybar =
                maxcamber * (1 - 2 * camberpose + 2 * camberpose * x - x^2) /
                (1 - camberpose)^2
        end
    else
        ybar = 0.0
    end
    return ybar
end

"""
    naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

Compute x, y airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

**Arguments:**
 - `c::Float` : Maximum camber value (percent of chord)
 - `p::Float` : Position along chord (in 10ths of chord) where maximum camber lies
 - `t::Float` : Maximum thickness of airfoil in percent chord

**Keyword Arguments:**
 - `N::Int` : Total number of coordinates to use (should be odd)
 - `x::Array{Float}` : x coordinates (cosine spaced coordinates used by default)
 - `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
 - `split::Bool` : Flag wheter to split into upper and lower halves.
"""
function naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if x == nothing
        x = cosine_spacing(N)
    end

    #naca digits
    maxcamber = c / 100
    camberpose = p / 10
    maxthick = t / 100

    #initialize arrays
    yu = zeros(N) #upper y values
    yl = zeros(N) #lower y values

    if real(c) != c || real(p) != p || real(t) != t
        yu = complex(yu)
        yl = complex(yl)
    end

    #--Calculate y-values--#
    for i in 1:N
        #thickness distribution
        T = thickness(x[i], maxthick; blunt_TE=blunt_TE)

        #camber distribution
        ybar = camber(x[i], maxcamber, camberpose)

        #y-positions at chordwise stations
        yl[i] = ybar - T / 2
        yu[i] = ybar + T / 2
    end

    if split
        return reverse(x), x, reverse(yl), yu
    else
        return [reverse(x); x[2:end]], [reverse(yl); yu[2:end]]
    end
end
