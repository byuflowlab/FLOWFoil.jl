"""
    whole_cosine_spacing(N::Integer=160)

Returns cosine spaced x coordinates from 1 to 0 back to 1.

# Arguments
- `N::Integer` : Total number of points is N+1.

# Returns
- `x::AbstractArray{Float}` : cosine spaced x-coordinates, starting at 1.0, going to 0.0, then back to 1.0.
"""
function whole_cosine_spacing(N::Integer=160)
    return [0.5 * (cos(2.0 * pi / (N) * (i - 1)) + 1) for i in 1:(N + 1)]
end

"""
    split_cosine_spacing(N::Integer=80)

Returns cosine spaced x coordinates from 0 to 1.

# Arguments
- `N::Integer` : Number of points.

# Returns
- `x::AbstractArray{Float}` : cosine spaced x-coordinates, starting at 0.0 ending at 1.0.
"""
function split_cosine_spacing(N::Integer=80)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
    repanel_airfoil(x, y; N=160)

Repanels airfoil coordinates using Akima splines with `N` coordinate points.

# Arguments
- `x::AbstractArray{Float}` : vector containing the x coordinates of the airfoil
- `y::AbstractArray{Float}` : vector containing the y coordinates of the airfoil

# Keyword Arguements
- `N::Int` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.

# Returns
- `repaneled_x::AbstractArray{Float}` : Repaneled, cosine spaced x corrdinates of the airfoil
- `repaneled_y::AbstractArray{Float}` : y coordinates of the repaneled airfoil obtained using an akima spline
"""
function repanel_airfoil(x, y; N=160)
    @assert length(x) == length(y) "x and y vectors must be the same length"

    #First normalize the airfoil to between 0 and 1
    normalize_coordinates!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N / 2)
    akimax = split_cosine_spacing(npoints)

    #now we split the top and bottom of the airfoil
    x1, x2, z1, z2 = split_upper_lower(x, y)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x1[1] > x1[end]
        x1 = reverse(x1)
        z1 = reverse(z1)
    end

    if x2[1] > x2[end]
        x2 = reverse(x2)
        z2 = reverse(z2)
    end

    #do the akima spline
    akimaz1 = FLOWMath.akima(x1, z1, akimax)
    akimaz2 = FLOWMath.akima(x2, z2, akimax)

    #figure out which spline is on top
    if maximum(akimaz1) > maximum(akimaz2) #then akimaz1 is on top so I need to reverse akimaz2
        repaneled_y = [reverse(akimaz2); akimaz1[2:end]]

    else #otherwise akimaz2 is on top so I need to reverse akimaz1
        repaneled_y = [reverse(akimaz1); akimaz2[2:end]]
    end

    repaneled_x = [reverse(akimax); akimax[2:end]]

    return repaneled_x, repaneled_y
end

"""
    repanel_airfoil(coordinates; N=160)

Repanels airfoil coordinates using Akima splines with `N` coordinate points.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates

# Keyword Arguements
- `N::Int=160` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.

# Returns
- `repaneled_coordinates::Array{Float}` : new coordinate array.
"""
function repanel_airfoil(coordinates; N=160)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, zpane = repanel_airfoil(x, y; N=N)

    return [xpane zpane]
end

"""
    repanel_revolution(coordinates; N=160)

Repanels body of revolution coordinates using Akima splines with `N` coordinate points.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates

# Keyword Arguements
- `N::Int=160` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.

# Returns
- `repaneled_coordinates::Array{Float}` : new coordinate array.
"""
function repanel_revolution(coordinates; N=160, normalize=true)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, ypane = repanel_revolution(x, y; N=N, normalize=normalize)

    return [xpane ypane]
end

"""
    repanel_revolution(x, y; N=160)

Repanels body of revolution coordinates using Akima splines with `N` coordinate points.

# Arguments
- `x::AbstractArray{Float}` : vector containing the x coordinates of the airfoil
- `y::AbstractArray{Float}` : vector containing the y coordinates of the airfoil

# Keyword Arguements
- `N::Int` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.

# Returns
- `repaneled_x::AbstractArray{Float}` : Repaneled, cosine spaced x corrdinates of the airfoil
- `repaneled_y::AbstractArray{Float}` : y coordinates of the repaneled airfoil obtained using an akima spline
"""

function repanel_revolution(x, y; N=160, normalize=true)
    @assert length(x) == length(y) "X and Y vectors must be the same length"

    if normalize
        scale = 1.0
        le = 0.0
    else
        scale = maximum(x) - minimum(x)
        le = minimum(x)
    end
    #First normalize the airfoil to between 0 and 1
    normalize_airfoil!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N)
    akimax = cosine_spacing(npoints)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x[1] > x[end]
        x = reverse(x)
        y = reverse(y)
    end

    #do the akima spline
    akimay = FLOWMath.akima(x, y, akimax)

    return akimax * scale .+ le, akimay * scale
end

"""
    refine_trailing_edge(coordinates)

Adds points along the trailing edge of an airfoil.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates

# Returns
- `x::AbstractArray{Float}` : vector containing the x coordinates of an airfoil with points along the trailing edge.
- `y::AbstractArray{Float}` : vector containing the y coordinates of an airfoil with points along the trailing edge.
"""
function refine_trailing_edge(coordinates)
    x = @view(coordinates[:, 1])
    y = @view(coordinates[:, 2])

    # Find appropriate trailing edge spacing
    # Take average of spacing of points adjacent to trailing edge.
    teapproxspacing =
        (
            sqrt((x[end - 1] - x[end])^2 + (y[end - 1] - y[end])^2) +
            sqrt((x[2] - x[1])^2 + (y[2] - y[1])^2)
        ) / 2
    # Find size of trailing edge
    tesize = abs(y[end] - y[1])
    # See how manz points is appropriate along trailing edge
    N = round(Int, tesize / teapproxspacing)
    if N > 1
        x, y = refine_trailing_edge(x, y, N)
    end
    return x, y
end

"""
    refine_trailing_edge(x, y, N::Integer)

Adds points along the trailing edge of an airfoil

# Arguments
- `x::AbstractArray{Float}` : vector containing the x coordinates of an airfoil.
- `y::AbstractArray{Float}` : vecotr containing the y coordinate of an  airfoil.
- `N::Integer` : Number of data points to be added at the trailing edge.

# Returns
- `x::AbstractArray{Float}` : vector containing the x coordinates of an airfoil with points along the trailing edge.
- `y::AbstractArray{Float}` : vector containing the y coordinates of an airfoil with points along the trailing edge.
"""
function refine_trailing_edge(x, y, N::Integer)
    TF = promote_type(eltype(x), eltype(y))
    x = vcat(x[1:(end - 1)], ones(TF, N - 1))
    zte = range(y[end], y[1], N)
    y = vcat(y[1:(end - 1)], zte[1:(end - 1)])
    return x, y
end

# """
#     fliplevelfix(x,y,angle,plot=true)
# Combines fixaf() rotate() fixaf() zerochordz() flipx() to level, flip,
# and normalize an airfoil for a given angle
# """
# function fliplevelfix(x, y, angle::Real)
#     x, y = fixaf(x, y)
#     x, y = rotateaf(x, y, angle)
#     x, y = fixaf(x, y)
#     x, y = flipx(x, y)
#     x, y = zerochordz(x, y)
#     return x, y
# end

#"""
#    fixaf(x,y)
#Adjusts coordinates of airfoil to loop from bottom edge trailing edge to top
#trailing edge
#"""
#function fixaf(coordinates)
#    x = coordinates[:, 1]
#    y = coordinates[:, 2]

#    # Ensure proper scaling
#    normalize_coordinates!(coordinates)

#    # Check if only top or bottom specified.  Assume symmetric and reflect if true.
#    if ((x[1] == 0.0) && (x[end] == 1.0)) || ((x[1] == 1.0) && (x[end] == 0.0))
#        idx = findall(y .!= 0.0)
#        xszm = x[idx]
#        zszm = -y[idx]
#        xszm = reverse(xszm; dims=1)
#        zszm = reverse(zszm; dims=1)
#        x = vcat(x, xszm)
#        y = vcat(y, zszm)
#    end

#    # Find trailing edge points
#    _, idxte = findmax(x)

#    #blunt trailing edge if trailing edge is at angle
#    if length(idxte) == 1 && y[idxte] != 0.0
#        blunt = true

#        # check which direction is more vertical
#        if idxte == length(x)
#            dotp1idx = 1
#        else
#            dotp1idx = idxte + 1
#        end

#        coordinates = [x y]
#        dotp1 =
#            abs(dot(coordinates[idxte, :] - coordinates[dotp1idx, :], [0.0, 1.0])) /
#            (norm(coordinates[idxte, :] - coordinates[dotp1idx, :]) * norm([0.0, 1.0]))

#        if idxte == 1
#            dotm1idx = length(x)
#        else
#            dotm1idx = idxte - 1
#        end

#        dotm1 =
#            abs(dot(coordinates[idxte, :] - coordinates[dotm1idx, :], [0.0, 1.0])) /
#            (norm(coordinates[idxte, :] - coordinates[dotm1idx, :]) * norm([0.0, 1.0]))

#        if dotp1 > dotm1
#            idxte = [idxte, dotp1idx]
#        else
#            idxte = [idxte, dotm1idx]
#        end

#    elseif length(idxte) > 1
#        blunt = true
#    else
#        blunt = false
#    end

#    idxstart = idxte[findmin(y[idxte])[2]]

#    # If trailing edge of pressure side is not first element, change accordingly
#    if (idxstart != 1)
#        x = vcat(x[idxstart:end], x[1:(idxstart - 1)])
#        y = vcat(y[idxstart:end], y[1:(idxstart - 1)])
#        idxte = 1
#    end

#    # Find leading edge
#    _, idxle = findmin(x)

#    # Check if not clockwise loop, flip if necessary
#    if (sum(y[1:idxle]) > sum(y[idxle:end]))
#        x = reverse(x; dims=1)
#        y = reverse(y; dims=1)
#        idxte = length(x)
#    end

#    # Check if blunt trailing edge and repeat te if not
#    if !blunt
#        x = vcat(x, x[idxte])
#        y = vcat(y, y[idxte])
#    end

#    return x, y
#end
