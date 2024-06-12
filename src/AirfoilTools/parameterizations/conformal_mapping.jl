
# TODO: update to use structs
# TODO: also make sure coordinate outputs are in standard format/order
##########################################################
##########################     ##########################
#####################     LOOK!    ######################
###########                                   ###########
#####     -----    TODO: YOU ARE HERE     -----     #####
###########                                   ###########
#####################     LOOK!    ######################
##########################     ##########################
#########################################################

######################################################################
#                                                                    #
#                           KARMAN TREFFTZ                           #
#                                                                    #
######################################################################

"""
"""
@kwdef struct KarmanTrefftz{Tb,Tr,Tw} <: AirfoilGeometry
    beta::Tb
    radius::Tr
    wedge_angle::Tw
end

"""
"""
function KarmanTrefftz(center, wedge_angle)
    return KarmanTrefftz(;
        radius=sqrt((1 - center[1])^2 + center[2]^2),
        beta=asin(center[2] / radius),
        wedge_angle=wedge_angle,
    )
end

"""
    karman_trefftz(paramters::KarmanTrefftz; N=360, normalize=true, split=false)

Karman-Trefftz airfoil parameterization based on angle beta, raidus, and wedge angle.

# Arguments:
 - `paramters::KarmanTrefftz` : KarmanTrefftz parameters

# Keyword Arguments:
 - `N::Int` : number of coordinates (for entire airfoil)
 - normalize::Bool` : Flag whether to normalize out put to unit chord and shift to have leading edge at zero.
 - split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
IF split == False:
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates
IF split == True:
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
"""
function karman_trefftz(p::KarmanTrefftz; N=360, normalize=true, split=false)
    return karman_trefftz(
        p.beta, p.radius, p.wedge_angle; N=N, normalize=normalize, split=split
    )
end

"""
    karman_trefftz(beta, radius, wedge; N=360, normalize=true, split=false)

Karman-Trefftz airfoil parameterization based on angle beta, raidus, and wedge angle.

# Arguments:
 - `beta::Float` : angle, in radians indicating center of circle relative to origin
 - `radius::Float` : radius of circle
 - `wedge::Float` : angle, in radians, of airfoil wedge angle

# Keyword Arguments:
 - `N::Int` : number of coordinates (for entire airfoil)
 - normalize::Bool` : Flag whether to normalize out put to unit chord and shift to have leading edge at zero.
 - split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
IF split == False:
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates
IF split == True:
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
"""
function karman_trefftz(beta, radius, wedge; N=360, normalize=true, split=false)

    #  convert wedge
    lambda = 2.0 - wedge / pi

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

    # separate x and z coordinates
    x = real(w)
    z = imag(w)

    if normalize
        normalize_coordinates!(x, z)
    end

    if split
        return split_upper_lower(x, z)
    else
        return x, z
    end
end

"""
    karman_trefftz(center, wedge; N=360, normalize=true, split=false)

Identical to `karman_trefftz(beta, radius, wedge)` but using center-based version.

# Arguments:
 - `center::Array{Float}` : [x z] location of circle center relative to origin
 - `wedge::Float` : angle, in radians, of airfoil wedge angle

# Keyword Arguments:
 - `N::Int` : number of coordinates (for entire airfoil)
 - `normalize::Bool` : Flag whether to normalize out put to unit chord and shift to have leading edge at zero.
 - `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
IF split == False:
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates
IF split == True:
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates
"""
function karman_trefftz(center, wedge; N=360, normalize=true, split=false)
    radius = sqrt((1 - center[1])^2 + center[2]^2) # radius

    beta = asin(center[2] / radius) # solve for beta

    return karman_trefftz(beta, radius, wedge; N=N, normalize=normalize, split=split)
end

######################################################################
#                                                                    #
#                             JOUKOWSKY                              #
#                                                                    #
######################################################################

"""
"""
@kwdef struct Joukowsky{Tc,Tr} <: AirfoilGeometry
    center::Tc
    radius::Tr
end

"""
    joukowsky(parameters::Joukowsky; N, fortest=false, normalize=true, split=false)

Joukowsky airfoil parameterization.

# Arguments:
 - `parameters::Joukowsky` : Joukowsky parameters

# Keyword Arguments:
 - `N::Int` : Number of coordinates to use
 - `fortest::Bool` : Flag to output non-coordinate paramters used in 'joukowsky_flow()'
 - `normalize::Bool` : Flag whether to normalize to unit chord and translate the leading edge to zero.
 - `split::Bool` : Flag wheter to split output into upper and lower surfaces.

# Returns:
IF split == False:
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates
IF split == True:
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates
"""
function joukowsky(p::Joukowsky; N=360, fortest=false, normalize=true, split=false)
    return joukowsky(
        p.center, p.radius; N=N, fortest=fortest, normalize=normalize, split=split
    )
end

"""
    joukowsky(center, radius; N, fortest=false, normalize=true, split=false)

Joukowsky airfoil parameterization.

# Arguments:
 - `center::Array{Float}` : [x z] location of center of circle relative to origin
 - `radius::Float` : radius of circle

# Keyword Arguments:
 - `N::Int` : Number of coordinates to use
 - `fortest::Bool` : Flag to output non-coordinate paramters used in 'joukowsky_flow()'
 - `normalize::Bool` : Flag whether to normalize to unit chord and translate the leading edge to zero.
 - `split::Bool` : Flag wheter to split output into upper and lower surfaces.

# Returns:
IF split == False:
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates
IF split == True:
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates
"""
function joukowsky(center, radius; N=360, fortest=false, normalize=true, split=false)
    beta = asin(center[2] / radius) # solve for beta

    theta = reverse(range(-beta; stop=2 * pi - beta, length=N))

    a = center[1] + radius * cos(beta)

    xi =
        (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta))) .+
        1.0 ./ (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta)))
    xi *= a

    x = real(xi)
    z = imag(xi)

    if normalize
        normalize_coordinates!(x, z)
    end

    if fortest
        return theta, beta, a, maximum(x) - minimum(x)
    else
        if split
            return split_upper_lower(x, z)
        else
            return x, z
        end
    end
end

"""
    joukowsky_flow(center, radius, alpha, vinf; N=360)

Calculate the analytic surface velocities and pressures as well as lift coefficient for a joukowsky airfoil.

# Arguments:
 - `center::Array{Float}` : [x z] location of circle center relative to origin
 - `radius::Float` : Radius of circle
 - `alpha::Float` : Angle of attack in degrees
 - `vinf::Float` : Freestream velocity

# Keyword Arguments:
 - `N::Int` : Number of coordinates to use

# Returns:
 - `vsurf::Array{Float}` : Magnitude of surface velocities at the nodes
 - `cpsurf::Array{Float}` : Surface pressures at the nodes
 - `cl::Float` : Lift coefficient
"""
function joukowsky_flow(center, radius, alpha, vinf; N=360)
    alpha_rad = alpha * pi / 180.0

    theta, beta, a, chord = joukowsky(center, radius; N=N, fortest=true)

    z = a * (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta)))

    vsurf =
        2.0 * vinf * (sin.(theta .- alpha_rad) .+ sin(alpha_rad .+ beta)) .* abs.(z) ./
        (abs.(z .- a^2 ./ z))

    cpsurf = 1.0 .- vsurf .^ 2.0 / vinf^2.0

    cl = 8.0 * pi * radius / chord * sin(alpha_rad + beta)

    return vsurf, cpsurf, cl
end
