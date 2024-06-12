"""
    z_from_parsec_coefficients(a, N::Int=80)

Calculate the x,z airfoil coordinates using the PARSEC polynomial.

# Arguments:
- `a::AbstractArray{Float}` : the PARSEC coefficients.
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

######################################################################
#                                                                    #
#                   MODIFIED PARSEC PARAMETERIZATION                 #
#                                                                    #
######################################################################

"""
# Fields:
- `leading_edge_radius::Float` : leading edge radius
- `maximum_thickness_xu::Float` : x-position of maximum thickness for upper side.
- `maximum_thickness_xl::Float` : x-position of maximum thickness for lower side.
- `maximum_thickness_zu::Float` : value of maximum thickness (from zero) for upper side.
- `maximum_thickness_zl::Float` : value of maximum thickness (from zero) for lower side.
- `curvature_u::Float` : curvature at point of maximum thickness on upper side.
- `curvature_l::Float` : curvature at point of maximum thickness on lower side.
- `trailing_edge_tangent_u::Float` : angle of surface at upper side trailing edge.
- `trailing_edge_tangent_l::Float` : angle of surface at lower side trailing edge.
- `trailing_edge_zu::Float=0.0` : z-position of upper side trailing edge.
- `trailing_edge_zl::Float=0.0` : z-position of lower side trailing edge.
"""
@kwdef struct ModifiedPARSEC{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11} <: AirfoilGeometry
    leading_edge_radius::T1
    maximum_thickness_xu::T2
    maximum_thickness_xl::T3
    maximum_thickness_zu::T4
    maximum_thickness_zl::T5
    curvature_u::T6
    curvature_l::T7
    trailing_edge_tangent_u::T8
    trailing_edge_tangent_l::T9
    trailing_edge_zu::T10 = 0.0
    trailing_edge_zl::T11 = 0.0
end

"""
    calculate_modified_parsec_coefficients(p, uppper_side)

Calculate the PARSEC coefficients using modified parameters (see parsec() docstring) for either the top or bottom curve.

# Arguments:
- `p::NamedTuple` : Named tuple of ModifiedPARSEC paramters including:
  - `leading_edge_radius::Float` : Leading edge radius
  - `maximum_thickness_x::Float` : chordwise position of maximum thickness
  - `maximum_thikcness_z::Float` : z-coordinate at maximum thickness
  - `curvature::Float` : second derivative of surface geometry at maximum thickness
  - `trailing_edge_tangent::Float` : trailing edge tangent angle, radians
  - `trailing_edge_z::Float` : z-position of trailing edge
- `side::Number` : +1 for upper side, -1 for lower side
"""
function calculate_modified_parsec_coefficients(p, side=1)
    TF = eltype(p)

    #---define b in Ax=b
    b = [
        sign(side) * sqrt(2 * p.leading_edge_radius)
        p.maximum_thickness_z
        0.0
        p.curvature
        tan(p.trailing_edge_tangent)
        p.trailing_edge_z
    ]

    #---define A in Ax=b
    #rename for convenience
    max_t_x = p.maximum_thickness_x

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
    r2 = [
        max_t_x^n2[1] max_t_x^n2[2] max_t_x^n2[3] max_t_x^n2[4] max_t_x^n2[5] max_t_x^n2[6]
    ]
    r3 = [
        c3[1]*max_t_x^n3[1] c3[2]*max_t_x^n3[2] c3[3]*max_t_x^n3[3] c3[4]*max_t_x^n3[4] c3[5]*max_t_x^n3[5] c3[6]*max_t_x^n3[6]
    ]

    r4 = [
        c4[1]*max_t_x^n4[1] c4[2]*max_t_x^n4[2] c4[3]*max_t_x^n4[3] c4[4]*max_t_x^n4[4] c4[5]*max_t_x^n4[5] c4[6]*max_t_x^n4[6]
    ]
    r5 = c6
    r6 = ones(1, 6)

    #assemble matrix
    A = [r1; r2; r3; r4; r5; r6]

    #solve for x in Ax=b (A\b)
    return A \ b
end

"""
    parsec(p::ModifiedPARSEC; N::Int=80, split=false)

Calculate the x,y airfoil coordinates for both top and bottom surfaces using modified PARSEC Parameterization method.

Use parsec() for standard PARSEC implementation.  This modified version employs direct values for trailing edge position and angles for each surface.

# Arguments:
- `p::ModifiedPARSEC` : ModifiedPARSEC paramters including:
  - `leading_edge_radius` : Leading edge radius
  - `maximum_thickness_xu` : chordwise position of maximum thickness of upper side
  - `maximum_thickness_xl` : chordwise position of maximum thickness of lower side
  - `maximum_thickness_zu` : z-coordinate at maximum thickness of upper side
  - `maximum_thickness_zl` : z-coordinate at maximum thickness of lower side
  - `curvature_u` : second derivative of surface geometry at maximum thickness of upper side
  - `curvature_l` : second derivative of surface geometry at maximum thickness of lower side
  - `trailing_edge_tangent_u` : trailing edge tangent angle of upper side
  - `trailing_edge_tangent_l` : trailing edge tangent angle of lower side
  - `trailing_edge_zu` : z-position of trailing edge of upper side
  - `trailing_edge_zl` : z-position of trailing edge of lower side

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function modified_parsec(p::ModifiedPARSEC; N::Integer=80, split=false)
    return modified_parsec(
        [
            p.leading_edge_radius
            p.maximum_thickness_xu
            p.maximum_thickness_xl
            p.maximum_thickness_zu
            p.maximum_thickness_zl
            p.curvature_u
            p.curvature_l
            p.trailing_edge_tangent_u
            p.trailing_edge_tangent_l
            p.trailing_edge_zu
            p.trailing_edge_zl
        ];
        N=N,
        split=split,
    )
end

"""
    parsec(p::AbstractArray{Float}; N::Int=80, split=false)

Calculate the x,y airfoil coordinates for both top and bottom surfaces using modified PARSEC Parameterization method.

Use parsec() for standard PARSEC implementation.  This modified version employs direct values for trailing edge position and angles for each surface.

# Arguments:
- `p::AbstractArray{Float}` : ModifiedPARSEC paramters including:
  - `leading_edge_radius` : Leading edge radius
  - `maximum_thickness_xu` : chordwise position of maximum thickness of upper side
  - `maximum_thickness_xl` : chordwise position of maximum thickness of lower side
  - `maximum_thickness_zu` : z-coordinate at maximum thickness of upper side
  - `maximum_thickness_zl` : z-coordinate at maximum thickness of lower side
  - `curvature_u` : second derivative of surface geometry at maximum thickness of upper side
  - `curvature_l` : second derivative of surface geometry at maximum thickness of lower side
  - `trailing_edge_tangent_u` : trailing edge tangent angle of upper side
  - `trailing_edge_tangent_l` : trailing edge tangent angle of lower side
  - `trailing_edge_zu` : z-position of trailing edge of upper side
  - `trailing_edge_zl` : z-position of trailing edge of lower side

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function modified_parsec(p; N::Integer=80, split=false)

    #--- Get x-values ---#
    x = split_cosine_spacing(N)

    #--- Upper Curve ---#
    au = calculate_modified_parsec_coefficients(
        (;
            leading_edge_radius=p[1],
            maximum_thickness_x=p[2],
            maximum_thickness_z=p[4],
            curvature=p[6],
            trailing_edge_tangent=p[8],
            trailing_edge_z=p[10],
        ),
        1,
    )
    zu = z_from_parsec_coefficients(au, N)

    #--- Lower Curve ---#
    al = calculate_modified_parsec_coefficients(
        (;
            leading_edge_radius=p[1],
            maximum_thickness_x=p[3],
            maximum_thickness_z=p[5],
            curvature=p[7],
            trailing_edge_tangent=p[9],
            trailing_edge_z=p[11],
        ),
        -1,
    )
    zl = z_from_parsec_coefficients(al, N)

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_modified_parsec(x,z)

Uses LsqFit to go from x-z coordinates to modified ModifiedPARSEC parameters.
"""
function determine_modified_parsec(x, z)
    function model(x, p)
        _, z = modified_parsec(p; N=ceil(Int, length(x) / 2))
        return z
    end

    #initial guess
    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

    fit = LsqFit.curve_fit(model, x, z, guess)

    return ModifiedPARSEC(;
        leading_edge_radius=fit.param[1],
        maximum_thickness_xu=fit.param[2],
        maximum_thickness_xl=fit.param[3],
        maximum_thickness_zu=fit.param[4],
        maximum_thickness_zl=fit.param[5],
        curvature_u=fit.param[6],
        curvature_l=fit.param[7],
        trailing_edge_tangent_u=fit.param[8],
        trailing_edge_tangent_l=fit.param[9],
        trailing_edge_zu=fit.param[10],
        trailing_edge_zl=fit.param[11],
    )
end

######################################################################
#                                                                    #
#                   STANDARD PARSEC PARAMETERIZATION                 #
#                                                                    #
######################################################################

"""
# Fields:
- `leading_edge_radius::Float` : leading edge radius
- `maximum_thickness_xu::Float` : x-position of maximum thickness for upper side.
- `maximum_thickness_xl::Float` : x-position of maximum thickness for lower side.
- `maximum_thickness_zu::Float` : value of maximum thickness (from zero) for upper side.
- `maximum_thickness_zl::Float` : value of maximum thickness (from zero) for lower side.
- `curvature_u::Float` : curvature at point of maximum thickness on upper side.
- `curvature_l::Float` : curvature at point of maximum thickness on lower side.
- `trailing_edge_angle::Float` : angle from chordline to horizontal at trailing edge.
- `boattail_angle::Float` : angle from chordline to upper/lower surfaces (half of wedge angle).
- `trailing_edge_gap::Float=0.0` : total gap distance between upper and lower surfaces at treailing edge.
- `trailing_edge_z::Float=0.0` : z-position of midpoint between upper and lower surfaces at trailing edge.
"""
@kwdef struct PARSEC{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11} <: AirfoilGeometry
    leading_edge_radius::T1
    maximum_thickness_xu::T2
    maximum_thickness_xl::T3
    maximum_thickness_zu::T4
    maximum_thickness_zl::T5
    curvature_u::T6
    curvature_l::T7
    trailing_edge_angle::T8
    boattail_angle::T9
    trailing_edge_gap::T10 = 0.0
    trailing_edge_z::T11 = 0.0
end

"""
    calculate_parsec_coefficients(p, side=1)

Calculate the PARSEC coefficients using standard parameters (see parsec() docstring) for either the top or bottom curve.

# Arguments:
- `p::NamedTuple` : NamedTuple of PARSEC standard paramters including:
  - [1]: `leading_edge_radius` : Leading edge radius
  - [2]: `X` : chordwise position of maximum thickness
  - [3]: `Z` : z-coordinate at maximum thickness
  - [4]: `curvature` : second derivative of surface geometry at maximum thickness
  - [5]: `trailing_edge_angle` : trailing edge angle
  - [6]: `boattail_angle` : boat-tail angle
  - [7]: `trailing_edge_gap` : z-distance between upper and lower surface trailing edge points
  - [8]: `trailing_edge_z` : z-position of center of trailing edge
"""
function calculate_parsec_coefficients(p, side=1)
    TF = eltype(p)

    #rename for convenience
    max_t_x = p.maximum_thickness_x

    # - define b in Ax=b - #
    b = [
        sign(side) * sqrt(2 * p.leading_edge_radius)
        p.trailing_edge_z + sign(side) * p.trailing_edge_gap / 2.0
        p.maximum_thickness_z
        tan(p.trailing_edge_angle + sign(side) * p.boattail_angle)
        0.0
        p.curvature
    ]

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
    r3 =
        [max_t_x^n2[1] max_t_x^n2[2] max_t_x^n2[3] max_t_x^n2[4] max_t_x^n2[5] max_t_x^n2[6]]
    r4 = c6
    r5 = [
        c3[1]*max_t_x^n3[1] c3[2]*max_t_x^n3[2] c3[3]*max_t_x^n3[3] c3[4]*max_t_x^n3[4] c3[5]*max_t_x^n3[5] c3[6]*max_t_x^n3[6]
    ]
    r6 = [
        c4[1]*max_t_x^n4[1] c4[2]*max_t_x^n4[2] c4[3]*max_t_x^n4[3] c4[4]*max_t_x^n4[4] c4[5]*max_t_x^n4[5] c4[6]*max_t_x^n4[6]
    ]

    # assemble matrix
    A = [r1; r2; r3; r4; r5; r6]

    # - solve for x in Ax=b (A\b) - #
    return A \ b
end

"""
    parsec(p::PARSEC; N::Integer=80, split=false)

Calculate the x,z airfoil coordinates for both top and bottom surfaces using standard PARSEC Parameterization method.

Use parsec() for modified PARSEC implementation.

# Arguments:
- `p::PARSEC` : PARSEC paramters including:
  - `leading_edge_radius` : Leading edge radius
  - `maximum_thickness_xu` : chordwise position of maximum thickness of upper side
  - `maximum_thickness_xl` : chordwise position of maximum thickness of lower side
  - `maximum_thickness_zu` : z-coordinate at maximum thickness of upper side
  - `maximum_thickness_zl` : z-coordinate at maximum thickness of lower side
  - `curvature_u` : second derivative of surface geometry at maximum thickness of upper side
  - `curvature_l` : second derivative of surface geometry at maximum thickness of lower side
  - `trailing_edge_angle` : trailing edge angle
  - `boattail_angle` : boat-tail angle
  - `trailing_edge_gap` : z-position of center of trailing edge
  - `trailing_edge_z` : z-distance between upper and lower surface trailing edge points

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function parsec(p::PARSEC; N::Integer=80, split=false)
    return parsec(
        [
            p.leading_edge_radius
            p.maximum_thickness_xu
            p.maximum_thickness_xl
            p.maximum_thickness_zu
            p.maximum_thickness_zl
            p.curvature_u
            p.curvature_l
            p.trailing_edge_angle
            p.boattail_angle
            p.trailing_edge_gap
            p.trailing_edge_z
        ];
        N=N,
        split=split,
    )
end

"""
    parsec(p::AbstractArray{Float}; N::Integer=80, split=false)

Calculate the x,z airfoil coordinates for both top and bottom surfaces using standard PARSEC Parameterization method.

Use parsec() for modified PARSEC implementation.

# Arguments:
- `p::AbstractArray{Float}` : PARSEC paramters including:
  - `leading_edge_radius` : Leading edge radius
  - `maximum_thickness_xu` : chordwise position of maximum thickness of upper side
  - `maximum_thickness_xl` : chordwise position of maximum thickness of lower side
  - `maximum_thickness_zu` : z-coordinate at maximum thickness of upper side
  - `maximum_thickness_zl` : z-coordinate at maximum thickness of lower side
  - `curvature_u` : second derivative of surface geometry at maximum thickness of upper side
  - `curvature_l` : second derivative of surface geometry at maximum thickness of lower side
  - `trailing_edge_angle` : trailing edge angle
  - `boattail_angle` : boat-tail angle
  - `trailing_edge_gap` : z-position of center of trailing edge
  - `trailing_edge_z` : z-distance between upper and lower surface trailing edge points

# Keyword Arguments:
- `N::Integer=80` : Number of x stations along chord
- `split::Bool` : Flag wheter to split into upper and lower halves.

# Returns:
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `z::AbstractArray{Float}` : Vector of z coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates from trailing edge to leading edge.
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates from leading edge to trailing edge.
"""
function parsec(p; N::Integer=80, split=false)

    #--- Get x-values ---#
    x = split_cosine_spacing(N)

    #--- Upper Curve ---#
    au = calculate_parsec_coefficients(
        (;
            leading_edge_radius=p[1],
            maximum_thickness_x=p[2],
            maximum_thickness_z=p[4],
            curvature=p[6],
            trailing_edge_angle=p[8],
            boattail_angle=p[9],
            trailing_edge_gap=p[10],
            trailing_edge_z=p[11],
        ),
        1,
    )
    zu = z_from_parsec_coefficients(au, N)

    #--- Lower Curve ---#
    al = calculate_parsec_coefficients(
        (;
            leading_edge_radius=p[1],
            maximum_thickness_x=p[3],
            maximum_thickness_z=p[5],
            curvature=p[7],
            trailing_edge_angle=p[8],
            boattail_angle=p[9],
            trailing_edge_gap=p[10],
            trailing_edge_z=p[11],
        ),
        -1,
    )
    zl = z_from_parsec_coefficients(al, N)

    if split
        return reverse(x), x, reverse(zl), zu
    else
        return [reverse(x); x[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    determine_parsec(x,z)

Uses LsqFit to go from x-z coordinates to standard PARSEC parameters.
"""
function determine_parsec(x, z)
    function model(x, p)
        _, z = parsec(p; N=ceil(Int, length(x) / 2))
        return z
    end

    #initial guess
    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

    fit = LsqFit.curve_fit(model, x, z, guess)

    return PARSEC(;
        leading_edge_radius=fit.param[1],
        maximum_thickness_xu=fit.param[2],
        maximum_thickness_xl=fit.param[3],
        maximum_thickness_zu=fit.param[4],
        maximum_thickness_zl=fit.param[5],
        curvature_u=fit.param[6],
        curvature_l=fit.param[7],
        trailing_edge_angle=fit.param[8],
        boattail_angle=fit.param[9],
        trailing_edge_gap=fit.param[10],
        trailing_edge_z=fit.param[11],
    )
end
