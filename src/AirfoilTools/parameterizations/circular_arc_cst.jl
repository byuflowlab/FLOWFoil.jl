"""
    CircularArcCST <: AirfoilGeometry

# Fields
- `maximum_camber::Float` : Value of maximimum camber in % chord
- `maximum_thickness_postition::Float` : Position of maximum thickness in % chord
- `maximum_thickness::Float` : Value of maximum thickness in % chord
- `A::Float=0.01` : Shape function parameter
- `k::Float=0.0` : Shape function parameter
- `z_te::Float=0.0` : Trailing edge thickness
"""
@kwdef struct CircularArcCST{TA,Tc,Tk,Tp,Tt,Tz} <: AirfoilGeometry
    maximum_camber::Tc
    maximum_thickness_postition::Tp
    maximum_thickness::Tt
    A::TA = 0.01
    k::Tk = 0.0
    z_te::Tz = 0.0
end

function circle_from_3pts(p1, p2, p3)
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    B = (x1^2 + y1^2) * (y3 - y2) + (x2^2 + y2^2) * (y1 - y3) + (x3^2 + y3^2) * (y2 - y1)
    C = (x1^2 + y1^2) * (x2 - x3) + (x2^2 + y2^2) * (x3 - x1) + (x3^2 + y3^2) * (x1 - x2)
    D = 2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2)

    xc = -B / D
    yc = -C / D
    R = sqrt((xc - x1)^2 + (yc - y1)^2)

    return xc, yc, R
end

function circular_arc_camber(maximum_camber, N)
    # Define three points on the arc
    A = [0.0, 0.0]
    B = [1.0, 0.0]
    M = [0.5, maximum_camber]

    xc, yc, R = circle_from_3pts(A, B, M)

    # Angles from center to start and end
    theta0 = atan(A[2] - yc, A[1] - xc)
    theta1 = atan(B[2] - yc, B[1] - xc)

    theta = linear_transform((0.0, 1.0), [theta0, theta1], split_cosine_spacing(N))
    x = xc .+ R .* cos.(theta)
    y = yc .+ R .* sin.(theta)

    #Ensure that the endpoints are where they should be
    x[1] = 0.0
    x[end] = 1.0
    y[1] = 0.0
    y[end] = 0.0

    return x, y, [xc, yc]
end

function class_fun(u)
    if iszero(u)
        return eltype(u)(0)
    else
        return sqrt(u) * (1.0 - u)
    end
end

function shape_fun(u, B, C; A=0.01, k=1.0)
    return A * cos(2.0 * pi * k * u) + B * u + C
end

function T(u, B, C; A=0.01, k=1.0, z_te=0.0)
    return shape_fun(u, B, C; A=A, k=k) * class_fun(u) + z_te * u
end

function residual!(r, y, x, p)
    maximum_thickness, maximum_thickness_postition = x
    (; A, k, z_te) = p
    B, C = y
    r[1] =
        T(maximum_thickness_postition, B, C; A=A, k=k, z_te=z_te) - maximum_thickness / 2.0
    r[2] =
        (
            0.5 * C + 1.5 * B * maximum_thickness_postition -
            1.5 * C * maximum_thickness_postition -
            2.5 * B * maximum_thickness_postition^2 +
            A *
            (0.5 - 1.5 * maximum_thickness_postition) *
            cos(2.0 * pi * k * maximum_thickness_postition) +
            A *
            k *
            maximum_thickness_postition *
            (-2.0 * pi + 2.0 * pi * maximum_thickness_postition) *
            sin(2.0 * pi * k * maximum_thickness_postition)
        ) / sqrt(maximum_thickness_postition) + z_te
    return r
end

function solve_BC(x, p)
    rwrap(r, y) = residual!(r, y, x, p)
    sol = nlsolve(rwrap, zeros(2); autodiff=:forward)
    return sol.zero
end

function get_normals(xc, yc, c)
    if isinf(c[2])
        nhat = [zeros(length(xc)) ones(length(xc))]
    else
        normals = [[x; y] .- c for (x, y) in zip(xc, yc)]
        nhat = vcat([n / norm(n) for n in normals]'...)
    end
    return nhat
end

function thicken_camber(xc, yc, tc, c)
    TF = promote_type(eltype(xc), eltype(yc), eltype(tc), eltype(c))
    upper_coordinates = zeros(TF, length(xc), 2)
    lower_coordinates = zeros(TF, length(xc), 2)

    normals = get_normals(xc, yc, c)

    for (i, (x, y, t, n)) in enumerate(zip(xc, yc, tc, eachrow(normals)))
        upper_coordinates[i, :] = [x; y] .+ t .* n
        lower_coordinates[i, :] = [x; y] .+ t .* (-n)
    end

    return upper_coordinates, lower_coordinates
end

"""
    circular_arc_cst(
        parameters::CircularArcCST; N=80, split=false, return_thickness_dist=false
    )

# Arguments
- `parameters::CircularArcCST` : Circular Arc CST Parameters

# Keyword Arguments
- `N::Int=161` : Total number of coordinates to use.  This values should be odd, but if not, the number of points returned will be N-1.
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `yl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
 - `yu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function circular_arc_cst(parameters::CircularArcCST; N=80, split=false)
    return circular_arc_cst(
        parameters.maximum_camber,
        parameters.maximum_thickness_postition,
        parameters.maximum_thickness,
        parameters.A,
        parameters.k,
        parameters.z_te;
        N=N,
        split=split,
    )
end

"""
    circular_arc_cst(
        maximum_camber,
        maximum_thickness_postition,
        maximum_thickness,
        A=0.01,
        k=0.0,
        z_te=0.0;
        N=80,
        split=false,
        return_thickness_dist=false,
    )

# Arguments
- `maximum_camber::Float` : Value of maximimum camber in % chord
- `maximum_thickness_postition::Float` : Position of maximum thickness in % chord
- `maximum_thickness::Float` : Value of maximum thickness in % chord
- `A::Float=0.01` : Shape function parameter
- `k::Float=0.0` : Shape function parameter
- `z_te::Float=0.0` : Trailing edge thickness

# Keyword Arguments
- `N::Int=161` : Total number of coordinates to use.  This values should be odd, but if not, the number of points returned will be N-1.
- `split::Bool=false` : Flag wheter to split into upper and lower halves.

# Returns
If `split` == false:
 - `x::AbstractArray{Float}` : Vector of x coordinates, clockwise from trailing edge.
 - `y::AbstractArray{Float}` : Vector of y coordinates, clockwise from trailing edge.
If `split` == true:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates from trailing edge to leading edge.
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates from leading edge to trailing edge.
 - `yl::AbstractArray{Float}` : Vector of lower half of y coordinates from trailing edge to leading edge.
 - `yu::AbstractArray{Float}` : Vector of upper half of y coordinates from leading edge to trailing edge.
"""
function circular_arc_cst(
    maximum_camber,
    maximum_thickness_postition,
    maximum_thickness,
    A=0.01,
    k=0.0,
    z_te=0.0;
    N=80,
    split=false,
)

    # get maximum_camber line
    xc, yc, c = circular_arc_camber(maximum_camber, N)

    # get coefficients for shape function
    B, C = iad.implicit(
        solve_BC,
        residual!,
        [maximum_thickness, maximum_thickness_postition],
        (; A=A, k=k, z_te=z_te),
    )

    # get thickness distribution
    tc = T.(xc, B, C; A=A, k=k, z_te=z_te)

    # place thickness along maximum_camber line
    upper_coordinates, lower_coordinates = thicken_camber(xc, yc, tc, c)

    if split
        return reverse(lower_coordinates[:, 1]),
        upper_coordinates[:, 1], reverse(lower_coordinates[:, 2]),
        upper_coordinates[:, 2]
    else
        return [reverse(lower_coordinates[:, 1]); upper_coordinates[2:end, 1]],
        [reverse(lower_coordinates[:, 2]); upper_coordinates[2:end, 2]]
    end
end
