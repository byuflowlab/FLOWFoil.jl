######################################################################
#                                                                    #
#                           KARMAN TREFFTZ                           #
#                                                                    #
######################################################################

"""
    KarmanTrefftz

# Fields
- `beta::Float` : angle, in radians indicating center of circle relative to origin
- `radius::Float` : radius of circle
- `wedge_angle::Float` : angle, in radians, of airfoil wedge angle
"""
@kwdef struct KarmanTrefftz{Tb,Tr,Tw} <: AirfoilGeometry
    beta::Tb
    radius::Tr
    wedge_angle::Tw
end

"""
    KarmanTrefftz(center, wedge_angle)

# Arguments
- `center::AbstractArray{Float}` : [x y] location of circle center relative to origin
- `wedge_angle::Float` : angle, in radians, of airfoil wedge angle

# Returns
- `parameters::KarmanTrefftz` : KarmanTrefftz parameters
"""
function KarmanTrefftz(center, wedge_angle)
    return KarmanTrefftz(;
        radius=sqrt((1 - center[1])^2 + center[2]^2),
        beta=asin(center[2] / radius),
        wedge_angle=wedge_angle,
    )
end

"""
    karman_trefftz(paramters::KarmanTrefftz; N=361, normalize=true, split=false)

Karman-Trefftz airfoil parameterization based on angle beta, raidus, and wedge angle.

# Arguments
- `paramters::KarmanTrefftz` : KarmanTrefftz parameters

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.
- `normalize::Bool=true` : Flag whether to normalize output to unit chord and shift to have leading edge at zero.
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false
- `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
- `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.

If `split` == true
- `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
- `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
- `yl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
- `yu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function karman_trefftz(p::KarmanTrefftz; N=361, normalize=true, split=false)
    return karman_trefftz(
        p.beta, p.radius, p.wedge_angle; N=N, normalize=normalize, split=split
    )
end

"""
    karman_trefftz(beta, radius, wedge_angle; N=361, normalize=true, split=false)

Karman-Trefftz airfoil parameterization based on angle beta, raidus, and wedge angle.

# Arguments
- `beta::Float` : angle, in radians indicating center of circle relative to origin
- `radius::Float` : radius of circle
- `wedge_angle::Float` : angle, in radians, of airfoil wedge angle

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.
- `normalize::Bool=true` : Flag whether to normalize output to unit chord and shift to have leading edge at zero.
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false
- `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
- `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.

If `split` == true
- `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
- `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
- `yl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
- `yu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function karman_trefftz(beta, radius, wedge_angle; N=361, normalize=true, split=false)

    #  convert wedge_angle
    lambda = 2.0 - wedge_angle / pi

    # if splitting, check that N is even
    if split && (N % 2 == 0.0)
        N += 1
        @warn("Number of desired points is not odd, adding one point for even split")
    end

    # initialize
    theta = range(-beta; stop=2 * pi - beta, length=N)

    # solve for geometry
    xi(theta) = 1.0 .+ radius * (exp(im * theta) - exp(-im * beta))
    inp(theta) = (1.0 .+ 1.0 / xi(theta))^(lambda)
    inm(theta) = (1.0 .- 1.0 / xi(theta))^(lambda)
    w = [
        lambda * (inp(theta[i]) + inm(theta[i])) / (inp(theta[i]) - inm(theta[i])) for
        i in 1:N
    ]

    # separate x and y coordinates
    x = real(w)
    y = imag(w)

    if normalize
        normalize_coordinates!(x, y)
    end

    if split
        return split_upper_lower(x, y)
    else
        return x, y
    end
end

"""
    karman_trefftz(center, wedge_angle; N=361, normalize=true, split=false)

Identical to `karman_trefftz(beta, radius, wedge_angle)` but using center-based version.

# Arguments
- `center::AbstractArray{Float}` : [x y] location of circle center relative to origin
- `wedge_angle::Float` : angle, in radians, of airfoil wedge angle

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.
- `normalize::Bool=true` : Flag whether to normalize out put to unit chord and shift to have leading edge at zero.
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
IF split == False
- `x::AbstractArray{Float}` : Array of x coordinates
- `y::AbstractArray{Float}` : Array of y coordinates

IF split == True
- `xu::AbstractArray{Float}` : Array of upper half of x coordinates
- `xl::AbstractArray{Float}` : Array of lower half of x coordinates
- `yu::AbstractArray{Float}` : Array of upper half of y coordinates
- `yl::AbstractArray{Float}` : Array of lower half of y coordinates
"""
function karman_trefftz(center, wedge_angle; N=361, normalize=true, split=false)
    radius = sqrt((1 - center[1])^2 + center[2]^2) # radius

    beta = asin(center[2] / radius) # solve for beta

    return karman_trefftz(beta, radius, wedge_angle; N=N, normalize=normalize, split=split)
end
