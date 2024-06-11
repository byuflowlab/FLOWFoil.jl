"""
    whole_cosine_spacing(N::Integer=160)

Returns cosine spaced x coordinates from 1 to 0 back to 1.

# Arguments:
- `N::Integer` : Total number of points is N+1.

# Returns:
- `x::Vector{Float}` : cosine spaced x-coordinates, starting at 1.0, going to 0.0, then back to 1.0.
"""
function whole_cosine_spacing(N::Integer=160)
    return [0.5 * (cos(2.0 * pi / (N) * (i - 1)) + 1) for i in 1:(N + 1)]
end

"""
    whole_cosine_spacing(N::Integer=160)

Returns cosine spaced x coordinates from 0 to 1.

# Arguments:
- `N::Integer` : Number of points.

# Returns:
- `x::Vector{Float}` : cosine spaced x-coordinates, starting at 0.0 ending at 1.0.
"""
function split_cosine_spacing(N::Integer=80)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
    repanel_airfoil(x,y;N)
    repanel_airfoil(xy;N)
Takes x and y coordinates of an airfoil  and uses a cosine spaced akima spline to fill in the gaps
**Arguments**
- `x::Vector{Float64}` : vector containing the x coordinates of the airfoil
- `y::Vector{Float64}` : vector containing the y components of the airfoil
- `xy::Array{Float64,2}` : Array of x and y coordinates with X in the first column and y in the 2nd
**Keyword Arguements**
- `N::Int` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.
**Returns**
- `xreturn::Vector{Float64}` : Repaneled, cosine spaced x corrdinates of the airfoil
- `yreturn::Vector{Float64}` : y coordinates of the repaneled airfoil obtained using an akima spline
- `xyreturn::Array{Float64}` : If the coordinates were input as an array, this will be returned with x in the 1st column and y in the 2nd.
"""
function repanel_airfoil(x, y; N=160)
    @assert length(x) == length(y) "X and Y vectors must be the same length"

    #First normalize the airfoil to between 0 and 1
    normalize_airfoil!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N / 2)
    akimax = cosine_spacing(npoints)

    #now we split the top and bottom of the airfoil
    x1, x2, y1, y2 = split_upper_lower(x, y)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x1[1] > x1[end]
        x1 = reverse(x1)
        y1 = reverse(y1)
    end

    if x2[1] > x2[end]
        x2 = reverse(x2)
        y2 = reverse(y2)
    end

    #do the akima spline
    akimay1 = FLOWMath.akima(x1, y1, akimax)
    akimay2 = FLOWMath.akima(x2, y2, akimax)

    #figure out which spline is on top
    if maximum(akimay1) > maximum(akimay2) #then akimay1 is on top so I need to reverse akimay2
        yreturn = [reverse(akimay2); akimay1[2:end]]

    else #otherwise akimay2 is on top so I need to reverse akimay1
        yreturn = [reverse(akimay1); akimay2[2:end]]
    end

    xreturn = [reverse(akimax); akimax[2:end]]

    return xreturn, yreturn
end

function repanel_airfoil(coordinates; N=160)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, ypane = repanel_airfoil(x, y; N=N)

    return [xpane ypane]
end

"""
    addtepoints(x::Array{<:Real,1},y::Array{<:Real,1})
Adds points along the trailing edge of an airfoil
"""
function addtepoints(x::Array{<:Real,1},y::Array{<:Real,1})
  # Find appropriate trailing edge spacing
  # Take average of spacing of points adjacent to trailing edge.
  teapproxspacing = (sqrt((x[end-1]-x[end])^2+(y[end-1]-y[end])^2)+
    sqrt((x[2]-x[1])^2+(y[2]-y[1])^2))/2
  # Find size of trailing edge
  tesize = abs(y[end]-y[1])
  # See how many points is appropriate along trailing edge
  N = round(Integer,tesize/teapproxspacing)
  x,y = addtepoints(x,y,N)
  return x,y
end

"""
    addtepoints(x::Array{<:Real,1},y::Array{<:Real,1},N::Integer)
Adds points along the trailing edge of an airfoil
"""
function addtepoints(x::Array{<:Real,1},y::Array{<:Real,1},N::Integer)
  x = vcat(x[1:(end-1)],ones(N-1))
  yte = linspace(y[end],y[1],N)
  y = vcat(y[1:(end-1)],yte[1:(end-1)])
  return x,y
end

"""
    fliplevelfix(x::Array{<:Real,1},y::Array{<:Real,1},angle,plot=true)
Combines fixaf() rotate() fixaf() zerochordy() flipx() to level, flip,
and normalize an airfoil for a given angle
"""
function fliplevelfix(x::Array{<:Real,1},y::Array{<:Real,1},angle::Real)
  x,y = fixaf(x,y)
  x,y = rotateaf(x,y,angle)
  x,y = fixaf(x,y)
  x,y = flipx(x,y)
  x,y = zerochordy(x,y)
  return x,y
end

"""
    fixaf(x,y)
Adjusts coordinates of airfoil to loop from bottom edge trailing edge to top
trailing edge
"""
function fixaf(x::Array{<:Real,1},y::Array{<:Real,1})
  # Ensure proper scaling
  minx = minimum(x)
  x = x-minx
  maxx = maximum(x)
  x = x/maxx
  y = y/maxx

  # Check if only top or bottom specified.  Assume symmetric and reflect if true.
  if ((x[1] == 0.0) && (x[end] == 1.0)) ||
    ((x[1] == 1.0) && (x[end] == 0.0))
    idx = find(y.!=0.0)
    xsym = x[idx]
    ysym = -y[idx]
    xsym = flipdim(xsym,1)
    ysym = flipdim(ysym,1)
    x = vcat(x,xsym)
    y = vcat(y,ysym)
  end

  # Find trailing edge points
  idxte = indmax(x)

  #blunt trailing edge if trailing edge is at angle
  if length(idxte) == 1 && y[idxte] !=0.0
    blunt = true
    # check which direction is more vertical
    coord = hcat(x,y)
    if idxte == length(x)
      dotp1idx=1
    else
      dotp1idx=idxte+1
    end
    dotp1 = abs(dot(coord[idxte,:]-coord[dotp1idx,:],[0.0,1.0]))/(norm(coord[idxte,:]-coord[dotp1idx,:])*norm([0.0,1.0]))
    if idxte == 1
      dotm1idx=length(x)
    else
      dotm1idx=idxte-1
    end
    dotm1 = abs(dot(coord[idxte,:]-coord[dotm1idx,:],[0.0,1.0]))/(norm(coord[idxte,:]-coord[dotm1idx,:])*norm([0.0,1.0]))
    if dotp1 > dotm1
      idxte = [idxte,dotp1idx]
    else
      idxte = [idxte,dotm1idx]
    end
  elseif length(idxte) > 1
    blunt = true
  else
    blunt = false
  end

  idxstart = idxte[indmin(y[idxte])]
  # If trailing edge of pressure side is not first element, change accordingly
  if (idxstart !=1)
    x = vcat(x[idxstart:end],x[1:(idxstart-1)])
    y = vcat(y[idxstart:end],y[1:(idxstart-1)])
    idxte = 1
  end
  # Find leading edge
  idxle = indmin(x)
  # Check if not clockwise loop, flip if necessary
  if (sum(y[1:idxle]) > sum(y[idxle:end]))
    x = flipdim(x,1)
    y = flipdim(y,1)
    idxte = length(x)
  end
  # Check if blunt trailing edge and repeat te if not
  if !blunt
    x = vcat(x,x[idxte])
    y = vcat(y,y[idxte])
  end
  return x,y
end


