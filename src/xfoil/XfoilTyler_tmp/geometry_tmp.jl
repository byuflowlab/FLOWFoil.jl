import Roots

struct Airfoil{TF}
    name::String
    x::Vector{TF}
    y::Vector{TF}
    gapTE::Bool         # if true a gap TE panel is used
end
Base.eltype(::Airfoil{TF}) where TF = TF

abstract type AbstractPanel end

struct Panel{TF} <: AbstractPanel
    x1::TF
    y1::TF
    x2::TF
    y2::TF
    xbar::TF
    ybar::TF
    length::TF
    stheta::TF
    ctheta::TF
end

struct WakePanel{TF} <: AbstractPanel
    x1::TF
    y1::TF
    x2::TF
    y2::TF
    xbar::TF
    ybar::TF
    length::TF
    stheta::TF
    ctheta::TF
    tdir1::AbstractVector
    tdir2::AbstractVector
end

struct TEgeometry{TF} <: AbstractPanel
    x1::TF
    y1::TF
    x2::TF
    y2::TF
    xbar::TF
    ybar::TF
    length::TF
    stheta::TF
    ctheta::TF
    tdp::TF
    tcp::TF
    hTE::TF # trailing edge gap - I believe this refers to the length of the gap normal to the 
                # streamline at the TE. So if the streamline is entirely in the horizontal direction, 
                # this would be the vertical length of the panel
    dtdx::TF # slope of TE thickness
end
Base.eltype(::Panel{TF}) where TF = TF

# convenience function to access fields within an array of structs
# function Base.getproperty(obj::AbstractVector{<:Panel}, sym::Symbol)
#     return getfield.(Ref(obj), sym)
# end

"""
    get_airfoil_from_file(file::AbstractString; name="")

For now, assumes that the airfoil file starts at the trailing edge and goes clockwise 
    around the airfoil (bottom surface first, then top surface).
"""
function get_airfoil_from_file(file::AbstractString; name="")

    # read in file coordinates
    data = readdlm(file, Float64; skipstart=1)

    # numpoints = length(data[:,1])   
    # ind_check = Int(ceil(numpoints / 4))

    # if data[1,1] == 0 # starts at zero --> reorder to start at 1
    #     ind = findfirst(x -> x==1.0, data[:,1])
    #     if data[ind_check,2] > data[end-ind_check,2] #starts along the top surface from LE
    #         x = vcat(data[ind:end,1],data[1:ind-1,1])
    #         y = vcat(data[ind:end,2],data[1:ind-1,2])
    #     else #starts along bottom surface from LE
    #         x = vcat(reverse(data[1:ind,1]),reverse(data[ind-1:1,1]))
    #         y = vcat(data[ind:end,2],data[1:ind,2]) #repeat TE point
    #     end
    # elseif data[1,1] == 1.0 #starts at TE   
    #     if data[ind_check,2] > data[end-ind_check,2] #going along top surface instead of bottom surface --> reverse entire order
    #         x = reverse(data[:,1])
    #         y = reverse(data[:,2])
    #     else #going in right direction
    #         x = data[:,1]
    #         y = data[:,2]
    #     end
    # end

    x = data[:,1]
    y = data[:,2]

    if x[end] === x[1] && y[end] === y[1]
        gapTE = true
    else
        gapTE = false
    end

    return Airfoil(name, x, y, gapTE)
end


function get_panel_nodes(airfoil, numpanels; cosine_spacing_flag=false)

    x = airfoil.x
    y = airfoil.y

    if !cosine_spacing_flag
        return x, y
    else
        @assert x[1] == 1.0 && x[end] == 1.0 "airfoil coordinates need to start and end at the trailing edge"
        numpoints = length(x)      
        ind_check = Int(ceil(numpoints / 4)) #1/4 chord
        @assert y[ind_check] < y[end-ind_check] "airfoil coordinates need to go clockwise (bottom surface first)"

        # split into upper and lower surfaces
        ind_0 = findfirst(z->z==minimum(x), x)
        x_lower = x[1:ind_0]
        y_lower = y[1:ind_0]
        x_upper = x[ind_0:end]
        y_upper = y[ind_0:end]

        x_cosine = cosine_spacing(eltype(x), ceil(Int64,numpanels/2+1))

        # akima spline of each surface with new xs
        y_akima_lower = FLOWMath.akima(reverse(x_lower), reverse(y_lower), x_cosine)
        y_akima_upper = FLOWMath.akima(x_upper, y_upper, x_cosine)

        # put together properly (not repeating LE, in correct order)
        xaf = vcat(reverse(x_cosine), x_cosine[2:end])
        yaf = vcat(reverse(y_akima_lower), y_akima_upper[2:end])
        
        return xaf, yaf
    end
end


"""
    cosine_spacing(TF, npoints)

Returns npoints-length vector of cosine spacing between 0 and 1.
"""
function cosine_spacing(TF::Type, npoints::Int64)

    x = Vector{TF}(undef, npoints)
    for i in 1:npoints
        x[i] = 0.5 * (1 - cos(pi * (i - 1) / (npoints - 1)))
    end

    return x
end

""" 
    get_panels(xaf, yaf)

    Define geometry of the airfoil, number of panels, coordinates for  each panel, angle for each panel, etc.
"""
function get_panels(xaf::AbstractVector, yaf::AbstractVector)

    if eltype(yaf) != Float64
        xaf = convert(Vector{eltype(yaf)}, xaf)
    end

    xbars = (xaf[1:end-1] .+ xaf[2:end])/2
    ybars = (yaf[1:end-1] .+ yaf[2:end])/2

    lengths = sqrt.((xaf[2:end] .- xaf[1:end-1]) .^ 2 .+ (yaf[2:end] .- yaf[1:end-1]) .^ 2)

    sthetas = (yaf[2:end] .- yaf[1:end-1]) ./ lengths
    cthetas = (xaf[2:end] .- xaf[1:end-1]) ./ lengths

    return Panel{eltype(xaf)}.(xaf[1:end-1], yaf[1:end-1], xaf[2:end], yaf[2:end], xbars, ybars, lengths, sthetas, cthetas)
end

function get_panels(airfoil::Airfoil, numpanels::Int64; cosine_spacing_flag=false)

    xaf, yaf = get_panel_nodes(airfoil, numpanels; cosine_spacing_flag)

    return get_panels(xaf, yaf)
end

function get_panels(;c=4, p=4, t=12, numpanels=160)

    xaf, yaf = naca4(c, p, t; N=numpanels+1)

    return get_panels(xaf, yaf)
end

function get_TE_panel(panels_airfoil)

    x1 = panels_airfoil[end].x2
    y1 = panels_airfoil[end].y2
    x2 = panels_airfoil[1].x1
    y2 = panels_airfoil[1].y1

    xbar = (x1 + x2) / 2
    ybar = (y1 + y2) / 2

    TE_length = sqrt((x2 - x1)^2 + (y2 - y1)^2)

    stheta_TE = (y2 - y1) / TE_length
    ctheta_TE = (x2 - x1) / TE_length

    t̂TE, t1, t2 = get_t̂TE(x2, y2, panels_airfoil[2].x1, panels_airfoil[2].y1, panels_airfoil[end].x1, panels_airfoil[end].y1, x1, y1)
    p̂TE = get_p̂TE(x2, y2, x1, y1)
    tdp = dot(t̂TE, p̂TE)
    tcp = cross_2D(t̂TE, p̂TE)

    s = [x1,y1] - [x2,y2]
    hTE = -s[1] * t̂TE[2] + s[2] * t̂TE[1]
    dtdx = t1[1] * t2[2] - t2[1] * t1[2]

    return TEgeometry(x1, y1, x2, y2, xbar, ybar, TE_length, stheta_TE, ctheta_TE, tdp, tcp, hTE, dtdx)
end

function get_t̂TE(x1, y1, x2, y2, xn, yn, xnp1, ynp1)

    vec1 = [x1, y1] - [x2, y2]
    vec1 /= norm(vec1)

    vecn = [xnp1, ynp1] - [xn, yn]
    vecn /= norm(vecn)

    t̂_TE = vec1 + vecn

    return t̂_TE / norm(t̂_TE), vec1, vecn
end

function get_p̂TE(x1, y1, xnp1, ynp1)

    vec = [xnp1, ynp1] .- [x1, y1]

    return vec ./ norm(vec)
end

function cross_2D(x,y)
    return x[1]*y[2] - x[2]*y[1]
end

"""
    function cosinespacing(N=180)

Calculate N cosine spaced points.

**Arguments:**
 - `N::Int` : Number of points

"""
function cosinespacing(N=80)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end


"""
inputs:
x: corresponding x locations
y: corresponding y locations
Vt: vector of tangential velocities around airfoil
panels: vector of panel objects - used to get panel lengths

assumes x, y is ordered starting from trailing edge going clockwise
around the airfoil back to the trailing edge 
(the convention we used for Hess-Smith)

also assumes that the tangential velocity from the panel method is
clockwise: positive on lower surface is upstream
(also Hess-Smith convention)
"""
function parseairfoil(x, y, Vt, panels, operatingparameters; verbose=false)
    
    if verbose; println("idx: $idx"); end
    
    idx = findfirst(z->z > 0, Vt)[1]
    idx_l = [idx-1, 1] #this lets me edit the entries if needed (the vector size is unchanged)
    idx_u = [idx, length(Vt)]
    lengths = [p.length for p in panels]

    # find stagnation point with akima interpolation
    s = vcat(0.0, cumsum(lengths))

    s_stagnation = ImplicitAD.implicit(solve_stagnation_state_implicit, residual_stagnation_implicit, vcat(s, Vt), nothing)
    
    # # separate into upper and lower surfaces
    sl = -(s[idx-1:-1:1] .- s_stagnation)
    su = s[idx:end] .- s_stagnation

    if operatingparameters.xft_l < 1
        operatingparameters.sft_l[1] = FLOWMath.linear(x[idx-1:-1:1], sl, operatingparameters.xft_l)
    else
        operatingparameters.sft_l[1] = s[end] #big enough it shouldn't ever use this
    end

    if operatingparameters.xft_u < 1
        operatingparameters.sft_u[1] = FLOWMath.linear(x[idx-1:-1:1], su, operatingparameters.xft_u)
    else
        operatingparameters.sft_u[1] = s[end] #big enough it shouldn't ever use this
    end

    return s, su, idx_u, sl, idx_l, s_stagnation
end

function solve_stagnation_state_implicit(x, p) #x needs to have all the possible dual variables

    npts = Int(length(x) / 2)
    s = x[1:npts]
    Vt = x[npts+1:end]

    spline = FLOWMath.Akima(s, Vt)

    s_stagnation_guess = FLOWMath.akima(Vt, s, 0.0)
    
    s_stagnation = Roots.find_zero(spline, s_stagnation_guess)

    return s_stagnation
end

function residual_stagnation_implicit(y, x, p)

    npts = Int(length(x) / 2)
    s = x[1:npts]
    Vt = x[npts+1:end]

    spline = FLOWMath.Akima(s, Vt)

    return spline(y)
end

"""
Assumes ue is all positive to begin with (from state vector).
"""
function split_s(s, ue, idx_u, idx_l, idx_w, swake0)

    Vt = vcat(reverse(-ue[idx_l[1]:-1:idx_l[end]]), ue[idx_u[1]:idx_u[end]])

    # spline = FLOWMath.Akima(s, Vt)
    # s_stagnation, info = FLOWMath.brent(spline, s[idx_l[1]], s[idx_u[1]])
    # s_stagnation = FLOWMath.akima(Vt, s, 0.0)
    
    s_stagnation = FLOWMath.linear(Vt, s, 0.0) #(in their initial s_stag calculation, mfoil solves for it, but in the update just uses a linear interpolation like this)
    su = s[idx_u[1]:idx_u[end]] .- s_stagnation
    sl = -(s[idx_l[1]:-1:idx_l[end]] .- s_stagnation)
    sw = swake0 .- s_stagnation

    return su, sl, sw, s_stagnation
end

# """
#     geometric_points(d1, L, N)

# Return a vector of N geometrically-spaced points starting at 0, 
# with the first interval `d1` and ending at `L`.
# Each interval increases by a constant ratio.
# """
# function geometric_points(d1, L, N)

#     # Solve for the geometric ratio r:
#     # sum = d1 * (1 - r^(N-1)) / (1 - r) = L
#     # Rearranged: d1 * (1 - r^(N-1)) / (1 - r) = L
#     # We'll solve for r numerically

#     f(r) = d1 * (1 - r^(N-1)) / (1 - r) - L

#     r, info = FLOWMath.brent(f, 1e-8, 10)

#     # Now build the points
#     points = vcat(0.0, cumsum(d1 * r .^ (0:N-2)))

#     return points
# end

function geometric_spaced_points(dx0, L, Np)
    # spaces Np points geometrically from [0,L], with dx0 as first interval
    # INPUTS
    #   dx0 : first interval length
    #   L   : total domain length
    #   Np  : number of points, including endpoints at 0,L
    # OUTPUTS
    #   x   : point locations (Vector{Float64} of length Np)

    N = Np - 1  # number of intervals

    d = L / dx0

    a = N * (N - 1) * (N - 2) / 6
    b = N * (N - 1) / 2
    c = N - d

    disc = max(b^2 - 4a*c, 0.0)
    r = 1 + (-b + sqrt(disc)) / (2a)

    # Newton iteration
    for _ in 1:10
        R   = r^N - 1 - d * (r - 1)
        R_r = N * r^(N - 1) - d
        dr  = -R / R_r
        if abs(dr) < 1e-6
            break
        end
        r += dr
    end

    # Build coordinate vector
    dx = dx0 .* r .^ (0:N-1)
    x = vcat(0.0, cumsum(dx))
# error()
    return x
end

function get_relative_geometry_panel_point(xi, yi, xj1, xj2, yj1, yj2, ctheta, stheta, length)

    xstar = (xi - xj1) * ctheta + (yi - yj1) * stheta
    ystar = -(xi - xj1) * stheta + (yi - yj1) * ctheta

    beta = get_beta_mfoil(1, 2, xi, yi, xj1, xj2, yj1, yj2)

    r1 = sqrt((xi - xj1)^2 + (yi - yj1)^2)
    r2 = sqrt((xi - xj2)^2 + (yi - yj2)^2)

    logr1 = iszero(r1) ? 0.0 : log(r1)
    logr2 = iszero(r2) ? 0.0 : log(r2)

    return xstar, ystar, r1, r2, logr1, logr2, beta
end

#### CST functions from NeualFoil.jl

"""
    bernstein(r, n, x)

Bernstein Basis Function: `binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)`
"""
function bernstein(r, n, x)
    return binomial(n, r) .* x .^ r .* (1 .- x) .^ (n .- r)
end

"""
    half_cst(coefficients, x, dz, leading_edge_weight; N1=0.5, N2=1.0)

Determine y-coordinates of one side of an airfoil give coeffiecients and x coordinates.

# Arguments
- `coefficients::Vector{Float}` : Kulfan parameters
- `x::Vector{Float}` : x-coordinates (front to back)
- `dz::Float` : Trailing edge gap
- `leading_edge_weight::Float` : Kulfan leading edge modification weight

# Keyword Arguments
- `N1::Float=0.5` : Class function parameter for leading edge
- `N2::Float=1.0` : Class function parameter for trailing edge

# Returns
- `y::Vector{Float}` : y-coordinates
"""
function half_cst(coefficients, x, dz, leading_edge_weight; N1=0.5, N2=1.0)
    nb = length(coefficients) - 1

    # Get class values
    C = @. x^N1 * (1.0 - x)^N2

    # Initialize shape functions
    S = similar(x) .= 0

    # Populate shape functions
    for (i, c) in enumerate(coefficients)
        S += c * bernstein(i - 1, nb, x)
    end

    # determine nominal y-values
    y = @. C * S + x * dz

    # Kulfan leading edge modification
    y .+= leading_edge_weight .* x .* (1.0 .- x) .^ (length(coefficients) + 0.5)

    return y
end

"""
    cst(x, p; N1=0.5, N2=1.0)

Determine y-coordinates of one side of an airfoil give coeffiecients and x coordinates.

# Arguments
- `x::Vector{Float}` : x-coordinates (concatenated top and bottom)
- `p::Vector{Float}` : parameters including Kulfan parameters, leading edge weight, and trailing edge gap.

# Keyword Arguments
- `N1::Float=0.5` : Class function parameter for leading edge
- `N2::Float=1.0` : Class function parameter for trailing edge

# Returns
- `y::Vector{Float}` : y-coordinates associated with the x-coordinates
"""
function cst(x, p; N1=0.5, N2=1.0)
    # Extract Parameters
    np = Int((length(p) - 2) / 2)
    pu = p[1:np]
    pl = p[(np + 1):(np * 2)]
    leading_edge_weight = p[end - 1]
    dz = p[end]

    # Split halves
    nx = Int(length(x) / 2)
    xu = x[1:nx]
    xl = x[(nx + 1):end]

    # Get y-values
    yu = half_cst(pu, xu, dz / 2.0, leading_edge_weight)
    yl = half_cst(pl, xl, -dz / 2.0, leading_edge_weight)

    # Return
    return [yu; yl]
end

function get_beta_mfoil(i, j, xi, yi, panel_j::AbstractPanel)

    return get_beta_mfoil(i, j, xi, yi, 
                            panel_j.x1, 
                            panel_j.x2, 
                            panel_j.y1, 
                            panel_j.y2)
end

function get_beta_mfoil(i, j, xi, yi, x1, x2, y1, y2)

    if i == j
        βij = π
    elseif i == j + 1 || (xi == x1 && yi == y1) || (xi == x2 && yi == y2)
        βij = 0.0
    else 
        num = (x1 - xi) * (y2 - yi) - (y1 - yi) * (x2 - xi)
        den = (x1 - xi) * (x2 - xi) + (y1 - yi) * (y2 - yi)

        βij = atan(num,den)
    end

    return βij
end

