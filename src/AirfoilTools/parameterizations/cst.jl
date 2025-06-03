"""
    CST

# Fields
- `upper_coefficients::AbstractArray{Float}` : Vector of coefficients defining the upper side.
- `lower_coefficients::AbstractArray{Float}` : Vector of coefficients defining the lower side.
- `trailing_edge_zu::Float=0.0` : z-position of the upper side trailing edge.
- `trailing_edge_zl::Float=0.0` : z-position of the lower side trailing edge.
- `N1::Float=0.5` : inherent parameter for round-nosed airfoils.
- `N2::Float=1.0` : inherent parameter for sharp trailing edge (with optional blunt trailing edge) airfoils.
"""
@kwdef struct CST{T1,T2,T3,T4,T5,T6}
    upper_coefficients::T1
    lower_coefficients::T2
    trailing_edge_zu::T5 = 0.0
    trailing_edge_zl::T6 = 0.0
    N1::T3 = 0.5
    N2::T4 = 1.0
end

"""
    bernstein(r, n, x)

Bernstein Basis Function: `binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)`
"""
function bernstein(r, n, x)
    return binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)
end

"""
    cst(
        parameters::CST;
        N::Integer=80,
        x=split_cosine_spacing(N),
        split=false,
    )

Obtain airfoil coordiantes (clockwise from trailing edge) from the class shape transformation (CST) parameterization.

# Arguments
- `parameters::CST` : CST parameters for airfoil.

# Keyword Arguments
- `N::Integer=80` : number of points to use for each side
- `x::AbstractArray{Float}=split_cosine_spacing(N)` : x-coordinates to use.
- `trailing_edge_zu::Float=0.0` : upper side trailing edge gap
- `trailing_edge_zl::Float=0.0` : lower side trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2
- `split::Bool=false` : if true, returns upper and lower coordinates separately as xl, xu, zl, zu rather than just x, z.

# Returns
If `split` == false
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.

 If `split` == true
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function cst(p::CST; N::Integer=80, x=split_cosine_spacing(N), split=false)
    return cst(
        p.upper_coefficients,
        p.lower_coefficients;
        trailing_edge_zu=p.trailing_edge_zu,
        trailing_edge_zl=p.trailing_edge_zl,
        N1=p.N1,
        N2=p.N2,
        N=N,
        x=x,
        split=split,
    )
end

"""
    cst(
        upper_coefficients,
        lower_coefficients;
        N::Integer=80,
        x=split_cosine_spacing(N),
        trailing_edge_zu=0.0,
        trailing_edge_zl=0.0,
        N1=0.5,
        N2=1.0,
        split=false,
    )

Obtain airfoil coordiantes (clockwise from trailing edge) from the class shape transformation (CST) parameterization.

# Arguments
- `upper_coefficients::AbstractArray{Float}` : Vector of CST coefficients for upper side of airfoil.
- `lower_coefficients::AbstractArray{Float}` : Vector of CST coefficients for lower side of airfoil.

# Keyword Arguments
- `N::Integer=80` : number of points to use for each side
- `x::AbstractArray{Float}=split_cosine_spacing(N)` : x-coordinates to use.
- `trailing_edge_zu::Float=0.0` : upper side trailing edge gap
- `trailing_edge_zl::Float=0.0` : lower side trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2
- `split::Bool=false` : if true, returns upper and lower coordinates separately as xl, xu, zl, zu rather than just x, z.

# Returns
If `split` == false
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.

If `split` == true
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function cst(
    upper_coefficients,
    lower_coefficients;
    N::Integer=80,
    x=split_cosine_spacing(N),
    trailing_edge_zu=0.0,
    trailing_edge_zl=0.0,
    N1=0.5,
    N2=1.0,
    split=false,
)
    zu = half_cst(upper_coefficients, x, trailing_edge_zu, N1, N2)
    zl = half_cst(lower_coefficients, x, trailing_edge_zl, N1, N2)

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
- `coefficients::AbstractArray{Float}` : Vector of CST coefficients.

# Keyword Arguments:
- `dz::Float=0.0` : trailing edge gap
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2
"""
function half_cst(coefficients, x=cosine_spacing(N), dz=0.0, N1=0.5, N2=1.0)
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
        x,
        z;
        n_upper_coefficients::Integer=8,
        n_lower_coefficients::Integer=8,
        trailing_edge_zu=0.0,
        trailing_edge_zl=0.0,
        N1=0.5,
        N2=1.0,
    )

Determine best-fit CST parameters using a least squares solve.

# Arguments
- `x::AbstractArray{Float}` : vector of x-coordinates.
- `z::AbstractArray{Float}` : vector of z-coordinates.

# Keyword Arguments
- `n_upper_coefficients::Integer=8` : number of upper side coefficients to fit
- `n_lower_coefficients::Integer=8` : number of lower side coefficients to fit
- `trailing_edge_zu::Float=0.0` : z coordiante of the trailing edge upper surface
- `trailing_edge_zl::Float=0.0` : z coordinate of the trailing edge lower surface
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2

# Returns
- `paramters::CST` : CST paramters for airfoil.
"""
function determine_cst(
    x, z; n_upper_coefficients::Integer=8, n_lower_coefficients::Integer=8, N1=0.5, N2=1.0
)

    # - split coordinates - #

    # find leading edge index
    _, leid = findmin(x)

    # lower side
    xl = x[1:leid]
    zl = z[1:leid]

    # upper side
    xu = x[leid:end]
    zu = z[leid:end]

    return determine_cst(
        xl,
        xu,
        zl,
        zu;
        n_upper_coefficients=n_upper_coefficients,
        n_lower_coefficients=n_lower_coefficients,
        N1=N1,
        N2=N2,
    )
end

"""
    determine_cst(
        xl,
        xu,
        zl,
        zu;
        n_upper_coefficients::Integer=8,
        n_lower_coefficients::Integer=8,
        N1=0.5,
        N2=1.0,
    )

Determine best-fit CST parameters for upper and lower sides of airfoil using a least squares solve.

# Arguments
- `xl::AbstractArray{Float}` : vector of lower side x-coordinates.
- `xu::AbstractArray{Float}` : vector of upper side x-coordinates.
- `zl::AbstractArray{Float}` : vector of lower side z-coordinates.
- `zu::AbstractArray{Float}` : vector of upper side z-coordinates.

# Keyword Arguments
- `n_upper_coefficients::Integer=8` : number of upper side coefficients to fit
- `n_lower_coefficients::Integer=8` : number of lower side coefficients to fit
- `trailing_edge_zu::Float=0.0` : z coordiante of the trailing edge upper surface
- `trailing_edge_zl::Float=0.0` : z coordinate of the trailing edge lower surface
- `N1::Float=0.5` : Class shape parameter 1
- `N2::Float=1.0` : Class shape parameter 2

# Returns
- `parameters::CST` : CST paramters for airfoil.
"""
function determine_cst(
    xl,
    xu,
    zl,
    zu;
    n_upper_coefficients::Integer=8,
    n_lower_coefficients::Integer=8,
    trailing_edge_zu=0.0,
    trailing_edge_zl=0.0,
    N1=0.5,
    N2=1.0,
)

    # # models to fit
    # cstU(x, p) = half_cst(p[1:(end - 1)], xu, p[end], N1, N2)
    # cstL(x, p) = half_cst(p[1:(end - 1)], xl, p[end], N1, N2)

    # # initial guesses
    # upper_coefficients = ones(n_upper_coefficients)
    # lower_coefficients = -ones(n_lower_coefficients)

    # # solve for coefficients
    # ufit = LsqFit.curve_fit(cstU, xu, zu, [upper_coefficients; 0.0])
    # lfit = LsqFit.curve_fit(cstL, xl, zl, [lower_coefficients; 0.0])

    # return CST(;
    #     upper_coefficients=ufit.param[1:(end - 1)],
    #     lower_coefficients=lfit.param[1:(end - 1)],
    #     N1=N1,
    #     N2=N2,
    #     trailing_edge_zu=ufit.param[end],
    #     trailing_edge_zl=lfit.param[end],
    # )

    upper_coeffs = determine_half_cst(
        xu,
        zu;
        n_coefficients=n_upper_coefficients,
        trailing_edge_z=trailing_edge_zu,
        N1=N1,
        N2=N2,
    )

    lower_coeffs = determine_half_cst(
        xl,
        zl;
        n_coefficients=n_lower_coefficients,
        trailing_edge_z=trailing_edge_zl,
        N1=N1,
        N2=N2,
    )

    return CST(;
        upper_coefficients=upper_coeffs[1:(end - 1)],
        lower_coefficients=lower_coeffs[1:(end - 1)],
        N1=N1,
        N2=N2,
        trailing_edge_zu=upper_coeffs[end],
        trailing_edge_zl=lower_coeffs[end],
    )
end
"""
    determine_half_cst(x, z; n_coefficients=8, trailing_edge_z=0.0, N1=0.5, N2=1.0)

Determine best-fit CST parameters for half of the airfoil using a least squares solve.

# Arguments
- `x::AbstractArray{Float}` : vector of x coordinates for airfoil
- `z::AbstractArray{Float}` : vecotr of z coordinates for airfoil

# Keyword Arguments
- `n_coefficients::Int=8` : number of coefficients to fit
- `trailing_edge_z::Float=0.0` : z coordinate of trailing edge
- `N1::Float=0.5` : Class Shape Parameter 1
- `N2::Float=1.0` : Class Shape Parameter 2

# Returns
- `parameters::CST` : CST parameters for airfoil
"""
function determine_half_cst(
    x, z; n_coefficients::Integer=8, trailing_edge_z=0.0, N1=0.5, N2=1.0
)

    # models to fit
    cst_model(x, p) = half_cst(p[1:(end - 1)], x, p[end], N1, N2)

    # initial guesses
    coefficients = ones(n_coefficients)
    coefficients .*= z[2] > z[1] ? 1.0 : -1.0

    # solve for coefficients
    fit = LsqFit.curve_fit(cst_model, x, z, [coefficients; 0.0])

    return fit.param
end

"""
    neuralfoil_half_cst(coefficients, x, dz, leading_edge_weight; N1=0.5, N2=1.0)

Compute the y-coordinates for one side (typically upper or lower) of an airfoil using the Class-Shape Transformation (CST) method with optional Kulfan leading edge shaping.

# Arguments
- `coefficients::AbstractArray{Float}`: Vector of shape coefficients (Bernstein polynomial weights).
- `x::AbstractArray{Float}`: Vector of x-coordinates (chordwise locations), typically ranging from 0 to 1.
- `dz::Float`: Linear thickness offset (typically used for camber line blending).
- `leading_edge_weight::Float`: Weight for the Kulfan-style leading edge shaping term.

# Keyword Arguments
- `N1::Float=0.5`: Class Shape Parameter 1
- `N2::Float=1.0`: Class Shape Parameter 2

# Returns
- `y::AbstractArray{Float}`: Vector of y-coordinates defining the airfoil surface (upper or lower).
"""
function neuralfoil_half_cst(coefficients, x, dz, leading_edge_weight; N1=0.5, N2=1.0)
    C = @. x^N1 * (1.0 - x)^N2

    nb = length(coefficients) - 1

    S = similar(x) .= 0

    for (i, c) in enumerate(coefficients)
        S += c * bernstein(i - 1, nb, x)
    end

    y = @. C * S + x * dz

    # Kulfan leading edge modification
    y .+= leading_edge_weight .* x .* (1.0 .- x) .^ (length(coefficients) + 0.5)

    return y
end

"""
    neuralfoil_cst(x, p; N1=0.5, N2=1.0)

Compute the y-coordinates for an airfoil using CST parameters.

# Arguments
- `x::AbstractArray{Float}`: Vector of x-coordinates
- `p::AbstractArray{Float}` : Vector of CST parameters (first half is for the upper surface and second half is for the lower surface)

# Keyword Arguments
- `N1::Float=0.5`: Class Shape Parameter 1
- `N2::Float=1.0`: Class Shape Parameter 2

# Returns
- `y::AbstractArray{Float}`: Vector of y-coordinates defining the airfoil with first half the upper surface and second half the lower.
"""
function neuralfoil_cst(x, p; N1=0.5, N2=1.0)
    np = Int((length(p) - 2) / 2)
    pu = p[1:np]
    pl = p[(np + 1):(np * 2)]
    leading_edge_weight = p[end - 1]
    dz = p[end]

    nx = Int(length(x) / 2)
    xu = x[1:nx]
    xl = x[(nx + 1):end]

    yu = neuralfoil_half_cst(pu, xu, dz / 2.0, leading_edge_weight)
    yl = neuralfoil_half_cst(pl, xl, -dz / 2.0, leading_edge_weight)

    return [yu; yl]
end

"""
    determine_neuralfoil_cst(coordinates; n_coefficients=8, N1=0.5, N2=1.0)

Fit CST (Class-Shape Transformation) parameters to a set of airfoil coordinates using nonlinear least squares optimization.

# Arguments
- `coordinates::Array{Float}`: Array of [x z] airfoil coordinates.

# Keyword Arguments
- `n_coefficients::Int=8`: Number of Bernstein polynomial coefficients to use for each surface (upper and lower).
- `N1::Float=0.5`: Class Shape Parameter 1
- `N2::Float=1.0`: Class Shape Parameter 2

# Returns
- `cst_upper::Vector{Float}`: CST shape coefficients for the upper surface.
- `cst_lower::Vector{Float}`: CST shape coefficients for the lower surface.
- `cst_LE::Float`: Leading-edge shape weight (Kulfan term).
- `cst_TE::Float`: Trailing edge thickness (`dz`).
"""
function determine_neuralfoil_cst(coordinates; n_coefficients=8, N1=0.5, N2=1.0)

    # Pre-process coordinates
    # Normalize
    normalize_coordinates!(coordinates)

    # Split
    xl, xu, yl, yu = split_upper_lower(coordinates; idx=Int((size(coordinates, 1) + 1) / 2))
    reverse!(xl)
    reverse!(yl)
    te_z = yu[end] - yl[end]

    fit = LsqFit.curve_fit(
        neuralfoil_cst,
        [xu; reverse(xl)],
        [yu; reverse(yl)],
        [0.1 * ones(n_coefficients); -0.1 * ones(n_coefficients); 0.1; te_z];
        autodiff=:forwarddiff,
    )

    nc = Int((length(fit.param) - 2) / 2)
    cst_upper = fit.param[1:nc]
    cst_lower = fit.param[(nc + 1):(end - 2)]
    cst_LE = fit.param[end - 1]
    cst_TE = fit.param[end]

    return cst_upper, cst_lower, cst_LE, cst_TE
end
