######################################################################
#                                                                    #
#                   MODIFIED PARSEC PARAMETERIZATION                 #
#                                                                    #
######################################################################

"""
"""
@kwdef struct PARSEC{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11} <: AirfoilGeometry
    rLE::T1
    xu::T2
    xl::T3
    zu::T4
    zl::T5
    Zxxup::T6
    Zxxlo::T7
    θTEup::T8
    θTElo::T9
    ZTEup::T10 = 0.0
    ZTElo::T11 = 0.0
end

# TODO: re-do parsec functions to use the type rather than arrays so that the various values are more precisely defined (rather than p[i], p.param)
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
    z_from_parsec_coefficients(a, N::Int=80)

Calculate the x,z airfoil coordinates using the PARSEC polynomial.

# Arguments:
- `a::Vector{Float}` : the PARSEC coefficients.
"""
function z_from_parsec_coefficients(a, N::Integer=80)
    x = split_cosine_spacing(N)
    z = similar(x) .= 0

    if real(a) != a
        x = complex(x)
        z = complex(z)
    end

    for i in 1:N
        for j in 1:6
            z[i] += a[j] * x[i]^(j - 0.5)
        end #for coeffs
    end #for x stations

    return z
end

"""
    calculate_parsec_coefficients(p, uppper_side)

Calculate the PARSEC coefficients using modified parameters (see parsec() docstring) for either the top or bottom curve.

# Arguments:
- `p::Vector{Float}` : Vector of PARSEC paramters including:
  - [1]: `rLE` : Leading edge radius
  - [2]: `X` : chordwise position of maximum thickness
  - [3]: `Z` : z-coordinate at maximum thickness
  - [4]: `Zxx` : second derivative of surface geometry at maximum thickness
  - [5]: `θTE` : trailing edge tangent angle
  - [6]: `ZTE` : z-position of trailing edge
- `side::Number` : +1 for upper side, -1 for lower side
"""
function calculate_parsec_coefficients(p, side=1)
    TF = eltype(p)

    #---define b in Ax=b
    b = [sign(side) * sqrt(2 * p[1]); p[3]; 0.0; p[4]; tan(p[5]); p[6]]

    #---define A in Ax=b
    #rename for convenience
    x = p[2]

    #calculate exponenets and coefficients
    n2 = zeros(TF, 6)
    n3 = zeros(TF, 6)
    n4 = zeros(TF, 6)
    c3 = zeros(TF, 6)
    c4 = zeros(TF, 6)
    for i in 1:6
        n2[i] = i - 0.5
        n3[i] = i - 1.5
        n4[i] = i - 2.5
        c3[i] = i - 0.5
        c4[i] = (i - 0.5) * (i - 1.5)
    end
    c6 = n2'

    #assemble matrix rows
    r1 = [1.0 0.0 0.0 0.0 0.0 0.0]
    r2 = [x^n2[1] x^n2[2] x^n2[3] x^n2[4] x^n2[5] x^n2[6]]
    r3 =
        [c3[1] * x^n3[1] c3[2] * x^n3[2] c3[3] * x^n3[3] c3[4] * x^n3[4] c3[5] * x^n3[5] c3[6] *
                                                                                         x^n3[6]]
    r4 =
        [c4[1] * x^n4[1] c4[2] * x^n4[2] c4[3] * x^n4[3] c4[4] * x^n4[4] c4[5] * x^n4[5] c4[6] *
                                                                                         x^n4[6]]
    r5 = c6
    r6 = ones(1, 6)

    #assemble matrix
    A = [r1; r2; r3; r4; r5; r6]

    #solve for x in Ax=b (A\b)
    return A \ b
end

"""
    parsec(p; N::Int=80, split=false)

Calculate the x,y airfoil coordinates for both top and bottom surfaces using modified PARSEC Parameterization method.

Use parsec_standard() for standard PARSEC implementation.  This modified version employs direct values for trailing edge position and angles for each surface.

# Arguments:
- `p::Vector{Float}` : Vector of PARSEC paramters including:
  - [1]: `rLE` : Leading edge radius
  - [2]: `xu` : chordwise position of maximum thickness of upper side
  - [3]: `xl` : chordwise position of maximum thickness of lower side
  - [4]: `zu` : z-coordinate at maximum thickness of upper side
  - [5]: `zl` : z-coordinate at maximum thickness of lower side
  - [6]: `Zxxup` : second derivative of surface geometry at maximum thickness of upper side
  - [7]: `Zxxlo` : second derivative of surface geometry at maximum thickness of lower side
  - [8]: `θTEup` : trailing edge tangent angle of upper side
  - [9]: `θTElo` : trailing edge tangent angle of lower side
  - [10]: `ZTEup` : z-position of trailing edge of upper side
  - [11]: `ZTElo` : z-position of trailing edge of lower side

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.
"""
function parsec(p; N::Integer=80, split=false)
    if length(p) == 9
        p = [p; 0.0; 0.0] #ZTEup = ZTElo = 0
    elseif length(p) == 12
        p = p[1:(end - 1)]
    elseif length(p) != 11
        error(
            "Incorrect number of parameters, must have one of the following: \n9: Sharp TE \n11: Full Parameter Set \n12: Full Set + AoA",
        )
    end
    #--- Get x-values ---#
    x = split_cosine_spacing(N)

    #--- Upper Curve ---#
    au = calculate_parsec_coefficients([p[1]; p[2]; p[4]; p[6]; p[8]; p[10]], 1)
    zu = z_from_parsec_coefficients(au, N)

    #--- Lower Curve ---#
    al = calculate_parsec_coefficients([p[1]; p[3]; p[5]; p[7]; p[9]; p[11]], -1)
    zl = z_from_parsec_coefficients(al, N)

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_parsec(x,z)

Uses LsqFit to go from x-z coordinates to modified PARSEC parameters.
"""
function determine_parsec(x, z)
    function model(x, p)
        _, z = parsec(p; N=ceil(Int, length(x) / 2))
        return z
    end

    #initial guess
    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

    fit = LsqFit.curve_fit(model, x, z, guess)
    p = fit.param'

    return p
end

######################################################################
#                                                                    #
#                   STANDARD PARSEC PARAMETERIZATION                 #
#                                                                    #
######################################################################

"""
    calculate_parsec_coefficients_standard(p, side=1)

Calculate the PARSEC coefficients using standard parameters (see parsec() docstring) for either the top or bottom curve.

# Arguments:
- `p::Vector{Float}` : Vector of PARSEC paramters including:
  - [1]: `rLE` : Leading edge radius
  - [2]: `X` : chordwise position of maximum thickness
  - [3]: `Z` : z-coordinate at maximum thickness
  - [4]: `Zxx` : second derivative of surface geometry at maximum thickness
  - [5]: `αTE` : trailing edge angle
  - [6]: `βTE` : boat-tail angle
  - [7]: `ΔZTE` : z-position of center of trailing edge
  - [8]: `ZTE` : z-distance between upper and lower surface trailing edge points
"""
function calculate_parsec_coefficients_standard(p, side=1)
    TF = eltype(p)

    #rename for convenience
    x = p[2]

    # - define b in Ax=b - #
    b = [sign(side) * sqrt(2 * p[1]); p[8] + p[7] / 2.0; p[3]; tan(p[5] - p[6]); 0.0; p[4]]

    # - define A in Ax=b - #

    # calculate exponenets and coefficients
    n2 = zeros(TF, 6)
    n3 = zeros(TF, 6)
    n4 = zeros(TF, 6)
    c3 = zeros(TF, 6)
    c4 = zeros(TF, 6)
    for i in 1:6
        n2[i] = i - 0.5
        n3[i] = i - 1.5
        n4[i] = i - 2.5
        c3[i] = i - 0.5
        c4[i] = (i - 0.5) * (i - 1.5)
    end
    c6 = n2'

    # assemble matrix rows
    r1 = [1.0 0.0 0.0 0.0 0.0 0.0]
    r2 = ones(1, 6)
    r3 = [x^n2[1] x^n2[2] x^n2[3] x^n2[4] x^n2[5] x^n2[6]]
    r4 = c6
    r5 = [
        c3[1]*x^n3[1] c3[2]*x^n3[2] c3[3]*x^n3[3] c3[4]*x^n3[4] c3[5]*x^n3[5] c3[6]*x^n3[6]
    ]
    r6 = [
        c4[1]*x^n4[1] c4[2]*x^n4[2] c4[3]*x^n4[3] c4[4]*x^n4[4] c4[5]*x^n4[5] c4[6]*x^n4[6]
    ]

    # assemble matrix
    A = [r1; r2; r3; r4; r5; r6]

    # - solve for x in Ax=b (A\b) - #
    return A \ b
end

"""
    parsec_standard(p; N::Integer=80, split=false)

Calculate the x,z airfoil coordinates for both top and bottom surfaces using standard PARSEC Parameterization method.

Use parsec() for modified PARSEC implementation.

# Arguments:
- `p::Vector{Float}` : Vector of PARSEC paramters including:
  - [1]: `rLE` : Leading edge radius
  - [2]: `xu` : chordwise position of maximum thickness of upper side
  - [3]: `xl` : chordwise position of maximum thickness of lower side
  - [4]: `zu` : z-coordinate at maximum thickness of upper side
  - [5]: `zl` : z-coordinate at maximum thickness of lower side
  - [6]: `Zxxup` : second derivative of surface geometry at maximum thickness of upper side
  - [7]: `Zxxlo` : second derivative of surface geometry at maximum thickness of lower side
  - [8]: `αTE` : trailing edge angle
  - [9]: `βTE` : boat-tail angle
  - [10]: `ΔZTE` : z-position of center of trailing edge (Optional)
  - [11]: `ZTE` : z-distance between upper and lower surface trailing edge points (Optional)

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.
"""
function parsec_standard(p; N::Integer=80, split=false)
    if length(p) == 9
        p = [p; 0.0; 0.0] #ZTE = ΔZTE = 0
    elseif length(p) != 11
        error(
            "Incorrect number of parameters, must have one of the following: \n9: Sharp TE \n11: Full Parameter Set",
        )
    end
    #--- Get x-values ---#
    x = split_cosine_spacing(N)

    #--- Upper Curve ---#
    au = calculate_parsec_coefficients_standard(
        [p[1]; p[2]; p[4]; p[6]; p[8]; p[9]; p[10]; p[11]], 1
    )
    zu = z_from_parsec_coefficients(au, N)

    #--- Lower Curve ---#
    al = calculate_parsec_coefficients_standard(
        [p[1]; p[3]; p[5]; p[7]; p[8]; p[9]; p[10]; p[11]], -1
    )
    zl = z_from_parsec_coefficients(al, N)

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_parsec_standard(x,z)

Uses LsqFit to go from x-z coordinates to standard PARSEC parameters.
"""
function determine_parsec_standard(x, z)
    function model(x, p)
        _, z = parsec_standard(p; N=ceil(Int, length(x) / 2))
        return z
    end

    #initial guess
    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

    fit = LsqFit.curve_fit(model, x, z, guess)
    p = fit.param'

    return p
end
