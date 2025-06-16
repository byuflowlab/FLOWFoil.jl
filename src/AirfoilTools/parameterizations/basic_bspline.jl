"""
    BasicBSpline

# Fields
- `leading_edge_radius::Float` : leading edge radius
- `trailing_edge_camber_angle::Float` : trailing edge camber angle (angle of chordline from horizontal at trailing edge).
- `wedge_angle::Float` : Wedge angle (angle between upper and lower surfaces at trailing edge).
- `trailing_edge_gap::Float=0.0` : distance between upper and lower surfaces at trailing edge. A value of zero indicates a sharp trailing edge.
- `third_ctrlpt_position::Float=1.0/3.0` : the position of the third control point.  This is an inherent value in the parameterization and if changed, the other parameters will not behave as they are defined here.
"""
@kwdef struct BasicBSpline{T1,T2,T3,T4,T5}
    leading_edge_radius::T1
    trailing_edge_camber_angle::T2
    wedge_angle::T3
    trailing_edge_gap::T4 = 0.0
    third_ctrlpt_position::T5 = 1.0 / 3.0
end

"""
    basic_bspline(parameters::BasicBSpline; N=160, split=false, return_nurbs=false)

Obtain airfoil coordinates from a B-Spline parameterization method.

# Arguments
- `parameters::BasicBSpline` : BasicBSpline parameters.

# Keyword Arguments
- `N::Integer=160` : number of points to use when defining the airfoil
- `split::Bool` : flag whether to output upper and lower coordinates separately
- `return_nurbs::Bool` : flag whether to output spline knots and control points as well

# Returns
if split=false
- `x::AbstractArray{Float}` : x-coordinates from lower TE clockwise to upper TE
- `z::AbstractArray{Float}` : z-coordinates from lower TE clockwise to upper TE

if split=true
- `xu::AbstractArray{Float}` : array of x-coordinates for the upper half of the airfoil (LE to TE)
- `zu::AbstractArray{Float}` : array of z-coordinates for the upper half of the airfoil (LE to TE)
- `xl::AbstractArray{Float}` : array of x-coordinates for the lower half of the airfoil (LE to TE)
- `zl::AbstractArray{Float}` : array of z-coordinates for the lower half of the airfoil (LE to TE)

if return_nurbs=true, also return:
- `NURBSu::NURBS.NURBScurve` : upper spline object
- `NURBSl::NURBS.NURBScurve` : lower spline object
"""
function basic_bspline(p::BasicBSpline; N=160, split=false, return_nurbs=false)
    return basic_bspline(
        p.leading_edge_radius,
        p.trailing_edge_camber_angle,
        p.wedge_angle;
        N=N,
        trailing_edge_gap=p.trailing_edge_gap,
        third_ctrlpt_position=p.third_ctrlpt_position,
        split=split,
        return_nurbs=return_nurbs,
    )
end

"""
    gbs(leading_edge_radius, trailing_edge_camber_angle, wedge_angle; perturbations=nothing, trailing_edge_gap=0, degree=3, third_ctrlpt_position=1/3, weights=nothing, split=false, return_nurbs=false)

Obtain airfoil coordinates from a B-Spline parameterization method.

# Arguments
- `leading_edge_radius::Float` : Leading Edge Radius
- `trailing_edge_camber_angle::Float` : Trailing Edge Camber Angle (degrees)
- `wedge_angle::Float` : Wedge Angle (degrees)

# Keyword Arguments
- `trailing_edge_gap::Float=0` : Trailing Edge Gap
- `third_ctrlpt_position::Float=1/3` : The x postion of the third control point.
- `split::Bool=false` : flag whether to output upper and lower coordinates separately
- `return_nurbs::Bool=false` : flag whether to output spline object

# Returns
if split=false
- `x::AbstractArray{Float}` : x-coordinates from lower TE clockwise to upper TE
- `z::AbstractArray{Float}` : z-coordinates from lower TE clockwise to upper TE

if split=true
- `xu::AbstractArray{Float}` : array of x-coordinates for the upper half of the airfoil (LE to TE)
- `zu::AbstractArray{Float}` : array of z-coordinates for the upper half of the airfoil (LE to TE)
- `xl::AbstractArray{Float}` : array of x-coordinates for the lower half of the airfoil (LE to TE)
- `zl::AbstractArray{Float}` : array of z-coordinates for the lower half of the airfoil (LE to TE)

if return_nurbs=true, also return
- `NURBSu::NURBS.NURBScurve` : upper spline object
- `NURBSl::NURBS.NURBScurve` : lower spline object
"""
function basic_bspline(
    leading_edge_radius,
    trailing_edge_camber_angle,
    wedge_angle;
    N=161,
    trailing_edge_gap=0,
    third_ctrlpt_position=1.0 / 3.0,
    split=false,
    return_nurbs=false,
)

    # get spline knots and control points
    unurbs, lnurbs = definespline(
        leading_edge_radius,
        trailing_edge_camber_angle,
        wedge_angle,
        trailing_edge_gap,
        third_ctrlpt_position,
    )

    #get coordinates from spline definition
    xu, zu = get_spline_coordinates(unurbs; N=floor(Int, N / 2))
    xl, zl = get_spline_coordinates(lnurbs; N=floor(Int, N / 2))

    #return what is asked for
    if return_nurbs && split
        return xu, zu, xl, zl, unurbs, lnurbs
    elseif return_nurbs && !split
        [reverse(xl); xu[2:end]], [reverse(zl); zu[2:end]], unurbs, lnurbs
    elseif !return_nurbs && split
        return xu, zu, xl, zl
    else
        return [reverse(xl); xu[2:end]], [reverse(zl); zu[2:end]]
    end
end

"""
    definespline(leading_edge_radius, trailing_edge_camber_angle, wedge_angle, trailing_edge_gap, third_ctrlpt_position)

Calculate the x and y location of the control points, and weight them according to
the weighting vector. Also, provide generic knot vectors.

# Arguments
- leading_edge_radius::Float : Leading edge radius
- trailing_edge_camber_angle::Float : trailing edge camber angle in degrees
- wedge_angle::Float : The trailing edge wedge angle in degrees
- trailing_edge_gap::Float : The trailing edge gap
- third_ctrlpt_position::Float : x position of the 2nd and 6th control points in degree 3.

# Returns
- `knots::AbstractArray{Float}`: The knot vector of the curve.
- `controlpoints::AbstractArray{Tuple{Float,3}}`: The control point vector, where each point is a tuple of (x, y, w).
"""
function definespline(
    leading_edge_radius,
    trailing_edge_camber_angle,
    wedge_angle,
    trailing_edge_gap,
    third_ctrlpt_position,
)
    TF = promote_type(
        eltype(third_ctrlpt_position),
        eltype(leading_edge_radius),
        eltype(trailing_edge_camber_angle),
        eltype(wedge_angle),
        eltype(trailing_edge_gap),
    )

    #Convert to radians
    boattailangle = pi * wedge_angle / 360 #boattail angle is half of wedge angle
    trailing_edge_camber_angle *= pi / 180

    #--Initialize knot vector
    knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    #calculate unweighted controlpoints
    Pl = [
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, -sqrt(2.0 * leading_edge_radius) / 3.0, 0.0),
        SVector(
            third_ctrlpt_position,
            (
                -trailing_edge_gap +
                (2.0 * tan(trailing_edge_camber_angle - boattailangle) / 3.0)
            ),
            0.0,
        ),
        SVector(1.0, -trailing_edge_gap / 2.0, 0.0),
    ]

    Pu = [
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, sqrt(2.0 * leading_edge_radius) / 3.0, 0.0),
        SVector(
            third_ctrlpt_position,
            (
                trailing_edge_gap +
                (2.0 * tan(trailing_edge_camber_angle + boattailangle) / 3.0)
            ),
            0.0,
        ),
        SVector(1.0, trailing_edge_gap / 2.0, 0.0),
    ]

    return NURBS.NURBScurve(NURBS.NURB(3, knots, ones(TF, size(Pu, 1))), Pu),
    NURBS.NURBScurve(NURBS.NURB(3, knots, ones(TF, size(Pl, 1))), Pl)
end

"""
    get_spline_coordinates(nurbs; N=80)

Compute and return plottable coordinates from a NURBS (Non-Uniform Rational B-Spline) curve.

# Arguments
- `nurbs::NURBS`: A NURBS curve object with fields:
    - `knots::AbstractVector{Float64}`: Knot vector
    - `controlPoints::AbstractMatrix{Float64}`: Weighted control points (each column is a point)
    - `degree::Int`: Degree of the curve

# Keyword Arguments
- `N::Int=80`: Number of segments (returns `N + 1` sampled points)

# Returns
- `x::Vector{Float64}`: x-coordinates of the sampled curve
- `z::Vector{Float64}`: z-coordinates of the sampled curve
"""
function get_spline_coordinates(nurbs; N=80)

    # create parametric point array
    u = collect(
        promote_type(eltype.(nurbs.controlPoints)...), range(0.0, 1.0; length=N + 1)
    )

    # get points
    points = nurbs(u)

    # separate points
    return getindex.(points, 1), getindex.(points, 2)
end

#"""
#    determine_basic_bspline(x,z)

#TODO: the output NURBS x-coordinates are not going to line up with the input x-coordinates.
#Need to figure out how to best fit the NURBS curve.
#"""
#function determine_basic_bspline(x, z)
#    function model(x, p)
#        _, z = basic_bspline(p[1], p[2], p[3]; N=length(x), trailing_edge_gap=p[4])
#        return z
#    end

#    #initial guess
#    guess = [0.001, 2.0, 14.0, 0.0]

#    fit = LsqFit.curve_fit(model, x, z, guess)

#    return BasicBSpline(;
#        leading_edge_radius=fit.param[1],
#        trailing_edge_camber_angle=fit.param[2],
#        wedge_angle=fit.param[3],
#        trailing_edge_gap=fit.param[4],
#        third_ctrlpt_position=1.0 / 3.0,
#    )
#end
