function cosinespacing(N::Int=80)
    x = zeros(N)
    for i in 1:N
        x[i] = 0.5 * (1 - cos(pi * (i - 1) / (N - 1)))
    end
    return x
end

"""
    thickness(x::Float64, maxthick::Number; bluntTE::Bool=false)

Compute thickness at a given chord-normalized x-position by NACA 4-series thickness equations.
"""
function thickness(x::Float64, maxthick::Number; bluntTE::Bool=false)
    c2 = 0.3537
    if bluntTE == true
        c2 = 0.3516
    end
    return 10 *
           maxthick *
           (0.2969 * sqrt(x) - 0.1260 * x - c2 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

"""
    camber(x::Float64, maxcamber::Number, camberpose::Number)

Compute camber at a given chord-normalized x-position by NACA 4-series camber equations.
"""
function camber(x::Float64, maxcamber::Number, camberpose::Number)
    if real(maxcamber) != 0.0 && real(camberpose) != 0.0
        if x <= real(camberpose)
            ybar = maxcamber * (2 * camberpose * x - x^2) / camberpose^2
        else
            ybar =
                maxcamber * (1 - 2 * camberpose + 2 * camberpose * x - x^2) /
                (1 - camberpose)^2
        end
    else
        ybar = 0.0
    end
    return ybar
end

"""
    naca4(c::Number=2.0, p::Number=4.0, t::Number=12.0, N::Int=100;bluntTE::Bool=false)

Compute x, y airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

x, y-upper, and y-lower coordinates are output. All have both leading and trailing edge points.

keyword argument: bluntTE=false (default) controls whether the blunt or sharp trailing edge formulation is used.
"""
function naca4(
    c::Number=2.0,
    p::Number=4.0,
    t::Number=12.0,
    N::Integer=100;
    x=cosinespacing(N),
    bluntTE::Bool=false,
)

    #naca digits
    maxcamber = c / 100
    camberpose = p / 10
    maxthick = t / 100

    #initialize arrays
    N = length(x)
    yu = zeros(N) #upper y values
    yl = zeros(N) #lower y values

    if real(c) != c || real(p) != p || real(t) != t
        yu = complex(yu)
        yl = complex(yl)
    end
    #--Calculate y-values--#
    for i in 1:N
        #thickness distribution
        T = thickness(x[i], maxthick; bluntTE=bluntTE)

        #camber distribution
        ybar = camber(x[i], maxcamber, camberpose)

        #y-positions at chordwise stations
        yl[i] = ybar - T / 2
        yu[i] = ybar + T / 2
    end

    return x, yu, yl
end
