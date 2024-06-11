"""
"""
@kwdef struct CST{Tu,Tl}
    upper_coefficients::Tu
    lower_coefficients::Tl
end

#TODO: switch to using struct
#########################################################
##########################     ##########################
#####################     LOOK!    ######################
###########                                   ###########
#####     -----    TODO: YOU ARE HERE     -----     #####
###########                                   ###########
#####################     LOOK!    ######################
##########################     ##########################
#########################################################

"""
    bernstein(r, n, x)

Bernstein Basis Function: `binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)`
"""
function bernstein(r, n, x)
    return binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)
end

"""
    cst(
        coefficients,
        N::Integer=80,
        x=split_cosine_spacing(N),
        dzu=0.0,
        dzl=0.0,
        N1=0.5,
        N2=1.0,
        split=false,
    )

Obtain airfoil coordiantes (clockwise from trailing edge) from the class shape transformation (CST) parameterization.

# Arguments:
- `upper_coefficients::Vector{Float}` : Vector of CST coefficients for upper side of airfoil.
- `lower_coefficients::Vector{Float}` : Vector of CST coefficients for lower side of airfoil.

# Keyword Arguments:
- `N::Integer=80` : number of points to use for each side
- `x::Vector{Float}=split_cosine_spacing(N)` : x-coordinates to use.
- `dzu::Float=0.0` : upper side trailing edge gap
- `dzl::Float=0.0` : lower side trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2
- `split::Bool=false` : if true, returns upper and lower coordinates separately as xl, xu, zl, zu rather than just x, z.

# Returns:
- `x::Vector{Float}` : vector of x-coordinates.
- `z::Vector{Float}` : vector of z-coordinates.
"""
function cst(
    upper_coefficients,
    lower_coefficients,
    N::Integer=80,
    x=split_cosine_spacing(N),
    dzu=0.0,
    dzl=0.0,
    N1=0.5,
    N2=1.0,
    split=false,
)
    zu = half_cst(upper_coefficients, x, dzu, N1, N2)
    zl = half_cst(lower_coefficients, x, dzl, N1, N2)

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    half_cst(coefficients, x=cosine_spacing(N), dz=0.0, N1=0.5, N2=1.0)

Define upper or lower side of airfoil using CST parameterization.

# Arguments:
- `coefficients::Vector{Float}` : Vector of CST coefficients.

# Keyword Arguments:
- `dz::Float=0.0` : trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2
"""
function half_cst(coefficients, x=cosine_spacing(N), dz=0.0, N1=0.5, N2=1.0)
    n = length(coefficients)

    C = @. x^N1 * (1.0 - x)^N2

    nb = length(coefficients) - 1

    S = similar(x) .= 0

    for (i, c) in enumerate(coefficients)
        S += c * bernstein(i - 1, nb, x)
    end

    return @. C * S + x * dz
end

"""
    determine_cst(
        xl,
        xu,
        zl,
        zu;
        n_upper_coefficients::Integer=8,
        n_lower_coefficients::Integer=8,
        dzu=0.0,
        dzl=0.0,
        N1=0.5,
        N2=1.0,
    )

Determine best-fit CST parameters for upper and lower sides of airfoil using a least squares solve.

# Arguments:
- `xl::Vector{Float}` : vector of lower side x-coordinates.
- `xu::Vector{Float}` : vector of upper side x-coordinates.
- `zl::Vector{Float}` : vector of lower side z-coordinates.
- `zu::Vector{Float}` : vector of upper side z-coordinates.

# Keyword Arguments:
- `n_upper_coefficients::Integer=8` : number of upper side coefficients to fit
- `n_lower_coefficients::Integer=8` : number of lower side coefficients to fit
- `dzu::Float=0.0` : upper side trailing edge gap
- `dzl::Float=0.0` : lower side trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2

# Returns:
- `upper_coefficients::Vector{Float}` : Vector of best-fit CST coefficients for upper side of airfoil.
- `lower_coefficients::Vector{Float}` : Vector of best-fit CST coefficients for lower side of airfoil.
"""
function determine_cst(
    xl,
    xu,
    zl,
    zu;
    n_upper_coefficients::Integer=8,
    n_lower_coefficients::Integer=8,
    dzu=0.0,
    dzl=0.0,
    N1=0.5,
    N2=1.0,
)

    # models to fit
    cstU(x, upper_coefficients) = half_cst(upper_coefficients, xu, dzu, N1, N2)
    cstL(x, lower_coefficients) = half_cst(lower_coefficients, xl, dzl, N1, N2)

    # initial guesses
    upper_coefficients = ones(n_coefficients)
    lower_coefficients = -ones(n_coefficients)

    # solve for coefficients
    ufit = LsqFit.curve_fit(cstU, xu, zu, upper_coefficients)
    lfit = LsqFit.curve_fit(cstL, xl, zl, lower_coefficients)

    return ufit.param, lfit.param
end
