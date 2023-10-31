#=
NACA 4-series Parametrization Functions
=#

"""
    thickness4(x, maxthick; blunt_TE=false)

Compute thickness at a given chord-normalized x-position by NACA 4-series thickness4 equations.

**Arguments:**
 - `x::Float` : x position along chordline
 - `maxthick::Float` : Maximum thickness value

**Keyword Arguments:**
 - `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
"""
function thickness4(x, maxthick; blunt_TE=false)

    # change c2 coefficient base on whether trailing edge is to be blunt or not.
    c2 = blunt_TE ? 0.3516 : 0.3537

    # return thickness4 value
    return 10.0 *
           maxthick *
           (0.2969 * sqrt(x) - 0.1260 * x - c2 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

"""
    camber4(x, maxcamber, camberpose)

Compute camber at a given chord-normalized x-position by NACA 4-series camber equations.

**Arguments:**
 - `x::Float` : x position along chordline
 - `maxcamber::Float64` : Maximum camber value
 - `camberpose::Float64` : Position of maximum camber
"""
function camber4(x, maxcamber, camberpose)
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
    naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

Compute x, y airfoil coordinates for N nodes, based on NACA 4-Series Parameterization.

**Arguments:**
 - `c::Float` : Maximum camber value (percent of chord)
 - `p::Float` : Position along chord (in 10ths of chord) where maximum camber4 lies
 - `t::Float` : Maximum thickness of airfoil in percent chord

**Keyword Arguments:**
 - `N::Int` : Total number of coordinates to use (should be odd)
 - `x::Array{Float}` : x coordinates (cosine spaced coordinates used by default)
 - `blunt_TE::Bool` : Flag whether trailing edge is blunt or not
 - `split::Bool` : Flag wheter to split into upper and lower halves.
"""
function naca4(c=2.0, p=4.0, t=12.0; N=161, x=nothing, blunt_TE=false, split=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = cosine_spacing(N)
    end

    #naca digits
    maxcamber = c / 100.0
    camberpose = p / 10.0
    maxthick = t / 100.0

    #initialize arrays
    TF = promote_type(typeof(c), eltype(x))
    yu = zeros(TF, N) #upper y values
    yl = zeros(TF, N) #lower y values

    #--Calculate y-values--#
    for i in 1:N
        #thickness4 distribution
        T = thickness4(x[i], maxthick; blunt_TE=blunt_TE)

        #camber4 distribution
        ybar = camber4(x[i], maxcamber, camberpose)

        #y-positions at chordwise stations
        yl[i] = ybar - T / 2
        yu[i] = ybar + T / 2
    end

    if split
        return reverse(x), x, reverse(yl), yu
    else
        return [reverse(x); x[2:end]], [reverse(yl); yu[2:end]]
    end
end

######################################################################
#                                                                    #
#               NACA 65-series Compressor Airfoil                    #
#                                                                    #
######################################################################
"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 1
"""
function thickness65c(x; method="scaled")
    if method == "derived"
        tr =
            [
                0.0
                0.772
                0.932
                1.169
                1.574
                2.177
                2.647
                3.040
                3.666
                4.143
                4.503
                4.760
                4.924
                4.996
                4.963
                4.812
                4.530
                4.146
                3.682
                3.156
                2.584
                1.987
                1.385
                0.810
                0.306
                0.0
            ] * 1e-2
        ler = 0.687
    elseif method == "scaled"
        tr =
            [
                0.0
                0.752
                0.890
                1.124
                1.571
                2.222
                2.709
                3.111
                3.746
                4.218
                4.570
                4.824
                4.982
                5.057
                5.029
                4.870
                4.570
                4.151
                3.627
                3.038
                2.451
                1.847
                1.251
                0.749
                0.354
                0.150
            ] * 1e-2
        ler = 0.666
    else
        @error "no method $method, please choose scaled or derived."
    end
    tx =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2

    return FLOWMath.akima(tx, tr, x)
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function camber65c(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    a1y =
        [
            0.0
            0.25
            0.35
            0.535
            0.93
            1.580
            2.120
            2.585
            3.365
            3.98
            4.475
            4.86
            5.15
            5.355
            5.475
            5.515
            5.475
            5.355
            5.15
            4.86
            4.475
            3.98
            3.365
            2.585
            1.58
            0.0
        ] * 1e-2

    a1fine = FLOWMath.akima(a1x, a1y, x)

    return a1fine * clo
end

"""
From NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds" Table 2
"""
function slopes65c(clo, x)
    a1x =
        [
            0.0
            0.5
            0.75
            1.25
            2.5
            5.0
            7.5
            10.0
            15.0
            20.0
            25.0
            30.0
            35.0
            40.0
            45.0
            50.0
            55.0
            60.0
            65.0
            70.0
            75.0
            80.0
            85.0
            90.0
            95.0
            100.0
        ] * 1e-2
    dydx = [
        0.4212
        0.38875
        0.3477
        0.29155
        0.2343
        0.19995
        0.17485
        0.13805
        0.11030
        0.08745
        0.06745
        0.04925
        0.03225
        0.01595
        0.0
        -0.01595
        -0.03225
        -0.04925
        -0.06745
        -0.08745
        -0.11030
        -0.13805
        -0.17485
        -0.23430
    ]

    dydxfine = FLOWMath.akima(a1x[2:(end - 1)], dydx, x)

    return dydx * clo
end

"""
Assumes x in non-dimensional range [0.0,1.0]


Description from NACA Research Memorandum L51G31: "Systematic Two-dimensional Cascade Tests of NACA 65-Series Compressor Blades at Low Speeds:"

The 65-series compressor blade family is formed by combining a basic thickness form with cambered mean lines.
The basic thickness form used is the NACA 65(216)-010 thickness form with the ordinates increased by 0.0015 times the chordwise stations to provide slightly, increased thickness toward the trailing edge.
In the scaled case, it was not derived for 10-percent thickness but was scaled down from the NACA 65,2-016 airfoil.
The scaling procedure gives the best results whep it is restricted to maximum thickness changes of a few percent.
The NACA 65-010 basic thickness has also been derived.
These thickness forms differ slightly but are considered to be interchangeable.

The basic mean line used is the a=1.0 mean line.
The amount of camber is for the design lift coefficient for the isolated airfoil with cl_o of 1.0.
Both ordinates and slopes are scaled directly to obtain other cambers.
Cambered blade sections are obtained by applying the thickness perpendicular to the mean line at stations laid out along the chord line.
In the designation the camber is given by the first number after the dash in tenths of cl_o.
For example, the NACA 65-810 and NACA 65-(12)10 blade sections are cambered for cl_o = 0.8 and cl_o = 1.2, respectively.
"""
function naca65c(clo; method="scaled", N=161, x=nothing, split=false)

    # get x coordinates
    N = Int(ceil(N / 2))
    if isnothing(x)
        x = cosine_spacing(N)
    end

    t = thickness65c(x; method=method)
    c = camber65c(clo, x)
    s = slopes65c(clo, x)

    #y-positions at chordwise stations
    yl = c .- t
    yu = c .+ t

    if split
        return reverse(x), x, reverse(yl), yu
    else
        return [reverse(x); x[2:end]], [reverse(yl); yu[2:end]]
    end
end
