"""
"""
@kwdef struct NACA4{Tc,Tp,Tt} <: AirfoilGeometry
    max_camber::Tc = 2.0
    max_campber_pos::Tp = 4.0
    max_thickness::Tt = 12.0
end

"""
    thickness4(x, maxthick; blunt_TE=false)

Compute thickness at a given chord-normalized x-position by NACA 4-series thickness4 equations.

# Arguments:
- `x::Float` : x position along chordlin, markersize=3, markershape=:squaree
- `maxthick::Float` : Maximum thickness value

# Keyword Arguments:
- `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
"""
function thickness4(x, maxthick; blunt_TE=false)

    # change c2 coefficient base on whether trailing edge is to be blunt or not.
    c2 = blunt_TE ? 0.3516 : 0.3537

    # return thickness4 value
    return 10.0 *
           maxthick *
           (0.2969 * sqrt(x) - 0.1260 * x - c2 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

"""
    camber4(x, maxcamber, camberpose)

Compute camber at a given chord-normalized x-position by NACA 4-series camber equations.

# Arguments:
- `x::Float` : x position along chordline
- `maxcamber::Float64` : Maximum camber value
- `camberpose::Float64` : Position of maximum camber
"""
function camber4(x, maxcamber, camberpose)
    if real(maxcamber) != 0.0 && real(camberpose) != 0.0
        if x <= real(camberpose)
            zbar = maxcamber * (2 * camberpose * x - x^2) / camberpose^2
        else
            zbar =
                maxcamber * (1 - 2 * camberpose + 2 * camberpose * x - x^2) /
                (1 - camberpose)^2
        end
    else
        zbar = 0.0
    end
    return zbar
end

"""
    naca4(parameters::NACA4; N=161, x=nothing, blunt_TE=false, split=false)
    naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

Compute x, z airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

# Arguments:
- `parameters::NACA4` : NACA 4-series parameters
- `c::Float` : Maximum camber value (percent of chord)
- `p::Float` : Position along chord (in 10ths of chord) where maximum camber4 lies
- `t::Float` : Maximum thickness of airfoil in percent chord

# Keyword Arguments:
- `N::Int` : Total number of coordinates to use (should be odd)
- `x::Array{Float}` : x coordinates (cosine spaced coordinates used by default)
- `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
- `split::Bool` : Flag wheter to split into upper and lower halves.
"""
function naca4(parameters::NACA4; N=161, x=nothing, blunt_TE=false, split=false)
    return naca4(
        parameters.max_camber,
        parameters.max_campber_pos,
        parameters.max_thickness;
        N=N,
        x=x,
        blunt_TE=blunt_TE,
        split=split,
    )
end

function naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = split_cosine_spacing(N)
    end

    #naca digits
    maxcamber = c / 100.0
    camberpose = p / 10.0
    maxthick = t / 100.0

    #initialize arrays
    TF = promote_type(typeof(c), eltype(x))
    zu = zeros(TF, N) #upper z values
    zl = zeros(TF, N) #lower z values

    #--Calculate z-values--#
    for i in 1:N
        #thickness4 distribution
        T = thickness4(x[i], maxthick; blunt_TE=blunt_TE)

        #camber4 distribution
        zbar = camber4(x[i], maxcamber, camberpose)

        #z-positions at chordwise stations
        zl[i] = zbar - T / 2
        zu[i] = zbar + T / 2
    end

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_naca4(x,z)

Calculate NACA 4-series parameters based on input x,z coordinates.
"""
function determine_naca4(x, z; blunt_TE=false)

    # model to fit
    function model(x, param)
        x, z = naca4(param[1], param[2], param[3]; x=x, blunt_TE=blunt_TE)
        return z
    end

    # solve for coefficients
    ufit = LsqFit.curve_fit(model, x, z, [2.0, 4.0, 12.0])

    # unpack
    return ufit.param[1], ufit.param[2], ufit.param[3]
end
