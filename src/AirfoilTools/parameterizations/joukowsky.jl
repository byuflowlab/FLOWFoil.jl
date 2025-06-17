"""
    Joukowsky

# Fields
- `center::AbstractArray{Float}` : [x y] location of center of circle relative to origin
- `radius::Float` : radius of circle
"""
@kwdef struct Joukowsky{Tc,Tr} <: AirfoilGeometry
    center::Tc
    radius::Tr
end

"""
    joukowsky(parameters::Joukowsky; N=361, fortest=false, normalize=true, split=false)

Joukowsky airfoil parameterization.

# Arguments
- `parameters::Joukowsky` : Joukowsky parameters

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.
- `fortest::Bool=false` : Flag to output non-coordinate paramters used in 'joukowsky_flow()'
- `normalize::Bool=true` : Flag whether to normalize to unit chord and translate the leading edge to zero.
- `split::Bool=false` : Flag wheter to split output into upper and lower surfaces.

# Returns
IF split == false
- `x::AbstractArray{Float}` : Array of x coordinates
- `y::AbstractArray{Float}` : Array of y coordinates

IF split == true
- `xu::AbstractArray{Float}` : Array of upper half of x coordinates
- `xl::AbstractArray{Float}` : Array of lower half of x coordinates
- `yu::AbstractArray{Float}` : Array of upper half of y coordinates
- `yl::AbstractArray{Float}` : Array of lower half of y coordinates
"""
function joukowsky(p::Joukowsky; N=361, fortest=false, normalize=true, split=false)
    return joukowsky(
        p.center, p.radius; N=N, fortest=fortest, normalize=normalize, split=split
    )
end

"""
    joukowsky(center, radius; N=361, fortest=false, normalize=true, split=false)

Joukowsky airfoil parameterization.

# Arguments
- `center::AbstractArray{Float}` : [x y] location of center of circle relative to origin
- `radius::Float` : radius of circle

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.
- `fortest::Bool=false` : Flag to output non-coordinate paramters used in 'joukowsky_flow()'
- `normalize::Bool=true` : Flag whether to normalize to unit chord and translate the leading edge to zero.
- `split::Bool=false` : Flag wheter to split output into upper and lower surfaces.

# Returns
IF split == false
- `x::AbstractArray{Float}` : Array of x coordinates
- `y::AbstractArray{Float}` : Array of y coordinates

IF split == true
- `xu::AbstractArray{Float}` : Array of upper half of x coordinates
- `xl::AbstractArray{Float}` : Array of lower half of x coordinates
- `yu::AbstractArray{Float}` : Array of upper half of y coordinates
- `yl::AbstractArray{Float}` : Array of lower half of y coordinates
"""
function joukowsky(center, radius; N=361, fortest=false, normalize=true, split=false)
    beta = asin(center[2] / radius) # solve for beta

    theta = reverse(range(-beta; stop=2 * pi - beta, length=N))

    a = center[1] + radius * cos(beta)

    xi =
        (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta))) .+
        1.0 ./ (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta)))
    xi *= a

    x = real(xi)
    y = imag(xi)

    if normalize
        normalize_coordinates!(x, y)
    end

    if fortest
        return theta, beta, a, maximum(x) - minimum(x)
    else
        if split
            return split_upper_lower(x, y)
        else
            return x, y
        end
    end
end

"""
    joukowsky_flow(center, radius, alpha; N=361)

Calculate the analytic surface velocities and pressures as well as lift coefficient for a joukowsky airfoil.

# Arguments
- `center::AbstractArray{Float}` : [x y] location of circle center relative to origin
- `radius::Float` : Radius of circle
- `alpha::Float` : Angle of attack in degrees

# Keyword Arguments
- `N::Int=361` : Total number of coordinates to use. Can be even or odd, but it is recommended to be odd for a clear leading edge point.

# Returns
- `vsurf::AbstractArray{Float}` : Magnitude of surface velocities at the nodes
- `cpsurf::AbstractArray{Float}` : Surface pressures at the nodes
- `cl::Float` : Lift coefficient
"""
function joukowsky_flow(center, radius, alpha; N=161)
    alpha_rad = alpha * pi / 180.0

    theta, beta, a, chord = joukowsky(center, radius; N=N, fortest=true)

    y = a * (1.0 .+ radius / a * (exp.(im * theta) .- exp(-im * beta)))

    vsurf =
        2.0 * (sin.(theta .- alpha_rad) .+ sin(alpha_rad .+ beta)) .* abs.(y) ./
        (abs.(y .- a^2 ./ y))

    cpsurf = 1.0 .- vsurf .^ 2.0

    cl = 8.0 * pi * radius / chord * sin(alpha_rad + beta)

    return vsurf, cpsurf, cl
end
