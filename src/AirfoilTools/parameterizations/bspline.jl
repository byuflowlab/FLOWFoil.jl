#=
BSpline parameterization from Rajnarayan using variable spline degree and purturbations.
Not currently working.
TODOS:
- swap to using GBS struct like other parameterizations
- Need to updated and get the perturbations working
- Need to swap to NURBS.jl rather than Splines.jl
- Need to wait for degree elevation to be added to NURBS.jl as well
- Will want to write a function to fit airfoil geometries using degree elevevation and perturbations as is done by Rajnarayan
=#

#"""
#"""
#@kwdef struct GBS{T1,T2,T3,T4,T5,T6,T7,T8}
#    leading_edge_radius::T1
#    trailing_edge_camber_angle::T2
#    wedge_angle::T3
#    perturbations::T4
#    trailing_edge_gap::T5 = 0.0
#    spline_degree::T6 = 3
#    third_ctrlpt_position::T7 = 1.0 / 3.0
#    weights::T8 = nothing
#end

#"""
#    gbs(parameters::GBS; N=160, split=false, debug=false)

#Obtain airfoil coordinates from a B-Spline parameterization method.

## Arguments:
#- `parameters::GBS` : GBS parameters.

## Keyword Arguments:
#- `N::Integer=160` : number of points to use when defining the airfoil
#- `split::Bool` : flag whether to output upper and lower coordinates separately
#- `debug::Bool` : flag whether to output spline knots and control points as well

## Returns:
#default:
#- `x::Array{Float}` : x-coordinates from lower TE clockwise to upper TE
#- `y::Array{Float}` : y-coordinates from lower TE clockwise to upper TE

#if split=true
#- `xu::Array{Float}` : array of x-coordinates for the upper half of the airfoil (LE to TE)
#- `yu::Array{Float}` : array of y-coordinates for the upper half of the airfoil (LE to TE)
#- `xl::Array{Float}` : array of x-coordinates for the lower half of the airfoil (LE to TE)
#- `yl::Array{Float}` : array of y-coordinates for the lower half of the airfoil (LE to TE)

#if debug keyword argument, also output
#- `uknots::Array{Float}` : Array of upper side spline knots
#- `ucontrolpoints::Array{Array{Float}}` : Array of upper side spline control points
#- `lknots::Array{Float}` : Array of lower side spline knots
#- `lcontrolpoints::Array{Array{Float}}` : Array of lower side spline control points
#"""
#function gbs(p::GBS; N=160, split=false, debug=false)
#    return gbs(
#        p.leading_edge_radius,
#        p.trailing_edge_camber_angle,
#        p.wedge_angle;
#        N=N,
#        perturbations=p.perturbations,
#        trailing_edge_gap=p.trailing_edge_gap,
#        degree=p.spline_degree,
#        third_ctrlpt_position=p.third_ctrlpt_position,
#        weights=p.weights,
#        split=split,
#        debug=debug,
#    )
#end

#"""
#    gbs(leading_edge_radius, trailing_edge_camber_angle, wedge_angle; perturbations=nothing, trailing_edge_gap=0, degree=3, third_ctrlpt_position=1/3, weights=nothing, split=false, debug=false)

#Obtain airfoil coordinates from a B-Spline parameterization method.

## Arguments:
#- `leading_edge_radius::Float` : Leading Edge Radius
#- `trailing_edge_camber_angle::Float` : Trailing Edge Camber Angle (degrees)
#- `wedge_angle::Float` : Wedge Angle (degrees)

## Keyword Arguments:
#- `purturbations::Array{Array{Float}}` : advanced control of airfoil shape. Each entry includes an upper and lower perturbation value, [pu; pl]
#- `trailing_edge_gap::Float` : Trailing Edge Gap
#- `degree::Int` : degree of NURBS curve used (for degrees higher than 3, control points are found with degree elevation.)
#- `third_ctrlpt_position::Float` : The x postion of the third control point.
#- `weights::Array{Array{Float}}` : a vector of weights for the control points. The length of the weights vector should be one more than the degree. Each entry has an upper-lower weight pair. [wu_i; wl_i]
#- `split::Bool` : flag whether to output upper and lower coordinates separately
#- `debug::Bool` : flag whether to output spline knots and control points as well

## Returns:
#default:
#- `x::Array{Float}` : x-coordinates from lower TE clockwise to upper TE
#- `y::Array{Float}` : y-coordinates from lower TE clockwise to upper TE

#if split=true
#- `xu::Array{Float}` : array of x-coordinates for the upper half of the airfoil (LE to TE)
#- `yu::Array{Float}` : array of y-coordinates for the upper half of the airfoil (LE to TE)
#- `xl::Array{Float}` : array of x-coordinates for the lower half of the airfoil (LE to TE)
#- `yl::Array{Float}` : array of y-coordinates for the lower half of the airfoil (LE to TE)

#if debug keyword argument, also output
#- `uknots::Array{Float}` : Array of upper side spline knots
#- `ucontrolpoints::Array{Array{Float}}` : Array of upper side spline control points
#- `lknots::Array{Float}` : Array of lower side spline knots
#- `lcontrolpoints::Array{Array{Float}}` : Array of lower side spline control points
#"""
#function gbs(
#    leading_edge_radius,
#    trailing_edge_camber_angle,
#    wedge_angle;
#    N=160,
#    perturbations=nothing,
#    trailing_edge_gap=0,
#    degree=3,
#    third_ctrlpt_position=1.0 / 3.0,
#    weights=nothing,
#    split=false,
#    debug=false,
#)

#    #do some checks
#    if weights != nothing
#        @assert length(weights) == degree + 1
#    end

#    # get spline knots and control points
#    unurbs, lnurbs = definespline(
#        degree,
#        third_ctrlpt_position,
#        leading_edge_radius,
#        trailing_edge_camber_angle,
#        wedge_angle,
#        trailing_edge_gap,
#        weights,
#    )

#    ##if using perturbations
#    #TODO: needs to be redone
#    #if perturbations != nothing

#    #    #check that perturbations have even length
#    #    @assert mod(length(perturbations), 2) == 0

#    #    #extract perturbations for each side
#    #    uperturbations = getindex.(perturbations, 1)
#    #    lperturbations = getindex.(perturbations, 2)

#    #    #update knots and control points based on perturbations
#    #    uknots, ucontrolpoints = perturb(
#    #        knots, ucontrolpoints, degree; deltaz=uperturbations
#    #    )
#    #    lknots, lcontrolpoints = perturb(
#    #        knots, lcontrolpoints, degree; deltaz=lperturbations
#    #    )
#    #else

#    #otherwise, knots don't change.
#    uknots = unurbs.knots
#    lknots = unurbs.knots
#    # end

#    #get coordinates from spline definition
#    xu, yu = get_spline_coordinates(unurbs; N=N)
#    xl, yl = get_spline_coordinates(lnurbs; N=N)

#    #return what is asked for
#    if debug && split
#        return xu, yu, xl, yl, unurbs, lnurbs
#    elseif debug && !split
#        [reverse(xl); xu[2:end]], [reverse(yl); yu[2:end]], unurbs, lnurbs
#        lcontrolpoints
#    elseif !debug && split
#        return xu, yu, xl, yl
#    else
#        return [reverse(xl); xu[2:end]], [reverse(yl); yu[2:end]]
#    end
#end

#"""
#    definespline()

#Calculate the x and y location of the control points, and weight them according to
#the weighting vector. Also, provide generic knot vectors.

## Inputs:
#- degree : degree of NURBS curve used
#- third_ctrlpt_position : x position of the 2nd and 6th control points in degree 3.
#- leading_edge_radius : Leading edge radius
#- wedge_angle : The trailing edge wedge angle
#- trailing_edge_gap : The trailing edge gap
#- weights : vector weights for the control points, leave empty for unity.

## Keyword Arguments:

## Outputs:
#- knots : the knot vector of the curve
#- controlpoints : control point vector (x,y,w)
#"""
#function definespline(
#    degree,
#    third_ctrlpt_position,
#    leading_edge_radius,
#    trailing_edge_camber_angle,
#    wedge_angle,
#    trailing_edge_gap,
#    weights,
#)
#    #Convert to radians
#    boattailangle = pi * wedge_angle / 360 #boattail angle is half of wedge angle
#    trailing_edge_camber_angle *= pi / 180

#    #check that degree is at least 3
#    # @assert degree >= 3
#    @assert degree == 3

#    #--Initialize knot vector
#    knots = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

#    #calculate unweighted controlpoints
#    Pl = [
#        [0.0, 0.0],
#        [0.0, -sqrt(2.0 * leading_edge_radius) / 3.0],
#        [
#            third_ctrlpt_position,
#            (
#                -trailing_edge_gap +
#                (2.0 * tan(trailing_edge_camber_angle - boattailangle) / 3.0)
#            ),
#        ],
#        [1.0, -trailing_edge_gap],
#    ]

#    Pu = [
#        [0.0, 0.0],
#        [0.0, sqrt(2.0 * leading_edge_radius) / 3.0],
#        [
#            third_ctrlpt_position,
#            (
#                trailing_edge_gap +
#                (2.0 * tan(trailing_edge_camber_angle + boattailangle) / 3.0)
#            ),
#        ],
#        [1.0, trailing_edge_gap],
#    ]

#    #if weights aren't given, go with ones as a default
#    if weights == nothing
#        wu = [1.0 for i in 1:4]
#        wl = [1.0 for i in 1:4]
#    else
#        wu = getindex(weights, 1)
#        wl = getindex(weights, 2)
#    end

#    #if degree is greater than 3, use degree elevation to update knots and control points.
#    # if degree == 3
#    #return splines
#    return Splines.NURBS(degree, knots, wu, Pu), Splines.NURBS(degree, knots, wl, Pl)

#    #elseif degree > 3
#    #    #raise degree
#    #    #TODO: this whole thing needs to be rewritten probably.

#    #    #Weight Control Points
#    #    uppercontrolpoints = [[0.0 for i in 1:3] for j in 1:4]
#    #    lowercontrolpoints = [[0.0 for i in 1:3] for j in 1:4]

#    #    for i in 1:4
#    #        uppercontrolpoints[i][1] = wu[i] * Pu[i, 1]
#    #        uppercontrolpoints[i][2] = wu[i] * Pu[i, 2]
#    #        uppercontrolpoints[i][3] = wu[i]
#    #    end
#    #    for i in 1:4
#    #        lowercontrolpoints[i][1] = wl[i] * Pl[i, 1]
#    #        lowercontrolpoints[i][2] = wl[i] * Pl[i, 2]
#    #        lowercontrolpoints[i][3] = wl[i]
#    #    end

#    #    _, raise_knots, raised_uppercontrolpoints = Splines.degreeelevatecurve(
#    #        length(uppercontrolpoints[:, 1]) - 1,
#    #        3,
#    #        tempknots,
#    #        uppercontrolpoints,
#    #        degree - 3,
#    #    )
#    #    _, raised_knots, raised_lowercontrolpoints = Splines.degreeelevatecurve(
#    #        length(lowercontrolpoints[:, 1]) - 1,
#    #        3,
#    #        tempknots,
#    #        lowercontrolpoints,
#    #        degree - 3,
#    #    )

#    #         return Splines.NURBS(
#    #             degree,
#    #             raised_knots,
#    #             [raised_uppercontrolpoints[i][end] for i in length(raised_uppercontrolpoints)],
#    #             [
#    #                 raised_uppercontrolpoints[i][1:(end - 1)] /
#    #                 raised_uppercontrolpoints[i][end] for
#    #                 i in 1:length(raised_uppercontrolpoints)
#    #             ],
#    #         )
#    #         return Splines.NURBS(
#    #             degree,
#    #             raised_knots,
#    #             [raised_lowercontrolpoints[i][end] for i in length(raised_lowercontrolpoints)],
#    #             [
#    #                 raised_lowercontrolpoints[i][1:(end - 1)] /
#    #                 raised_lowercontrolpoints[i][end] for
#    #                 i in 1:length(raised_lowercontrolpoints)
#    #             ],
#    #         )
#    #     else
#    #         @error "No functionality for degree <3"
#    #     end
#end

#"""
#    perturb()

#Insert knots around a certain location in an existing spline.

## Arguments:
#- `knots::Array{Float} : knot vector before insertion of new knots
#- `controlpoints::Array{Float,2}` : weighted control point vector
#- `degree::Int` : the degree of the curve
#- `deltaz::Array{Float}` : Array of perturbations

## Keyword Arguments:
#- `knotstoadd::Array{Float}` : user define knots to add.

## Returns:
#- `perturbed_knots::Array{Float}` : new knot vector with knots inserted
#- `perturbed_controlpoints::Array{Float,2}` : new control points vector (x,y,w)
#"""
#function perturb(knots, controlpoints, degree, deltaz; knotstoadd=nothing)

#    #get lengths
#    numcps = length(controlpoints[:, 1])

#    #Calculate new knots to add
#    if knotstoadd == nothing
#        knotstoadd = collect(0.0:(1.0 / (length(deltaz) + 1.0)):1.0)[2:(end - 1)]
#    end

#    #Use splines package to insert new knots.
#    refined_knots, refined_controlpoints = Splines.refineknotvectorcurve(
#        length(controlpoints[:, 1]) - 1,
#        degree,
#        knots,
#        controlpoints,
#        knotstoadd,
#        length(knotstoadd) - 1,
#    )

#    #find indices in new knot vector associated with new knots
#    newidx = findfirst.(isequal.(knotstoadd), (refined_knots,))

#    #Perturb the inserted knots in the y direction.
#    for i in 1:length(deltaz)
#        refined_controlpoints[newidx[i], 2] += deltaz[i]
#    end

#    return refined_knots, refined_controlpoints
#end

#"""
#        get_spline_coordinates()
#Get plotable coordinates from the weighted control points.

## Arguments:
#- `knots::Array{Float}` : knot vector
#- `controlpoints::Array{Float,2} `: Control Point matrix
#- `degree::Int` : degree of the NURBS curve

## Keyword Arguments:
#- `N::Int` : number of coordinates, or panels that will be generated

## Returns:
#- `x::Array{Float}` : x Airfoil coordinates
#- `y::Array{Float}` : y Airfoil coordinates
#"""
#function get_spline_coordinates(nurbs; N=160)

#    #create parametric point array
#    u = range(0.0, 1.0; length=N + 1)

#    #n = number of control points - 1
#    n = length(nurbs.ctrlpts) - 1

#    #initialize curve point array
#    Cw = [[0.0; 0.0] for i in 1:length(u)]

#    #loop through parametric points to get curve points
#    for i in 1:length(u)
#        Cw[i] = Splines.curvepoint(nurbs, u[i])
#    end

#    return getindex.(Cw, 1), getindex.(Cw, 2)
#end

#"""
#    determine_gbs(x,y)

#Uses LsqFit to go from xy to modified parsec without calculating each parameter.
#Depending on the initial guess and subsequent iterations, it may throw an error
#involving complex numbers. If so, alter the initial guess.
#"""
#function determine_bbs(x, y)
#    function model(x, p)
#        _, y = gbs(p[1], p[2], p[3]; N=ceil(Int, length(x) / 2))
#        return y
#    end

#    #initial guess
#    guess = [0.015, 0.0, 14.0]

#    fit = LsqFit.curve_fit(model, x, y, guess)

#    return GBS(;
#        leading_edge_radius=fit.param[1],
#        trailing_edge_camber_angle=fit.param[2],
#        wedge_angle=fit.param[3],
#        trailing_edge_gap=fit.param[4],
#        third_ctrlpt_position=1.0 / 3.0,
#        weights=nothing,
#    )
#end
