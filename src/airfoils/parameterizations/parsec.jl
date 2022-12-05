import NLsolve
import Dierckx

"""
    getxy(a::Array, N::Int=80)

Calculate the x,y airfoil coordinates using the PARSEC polynomial.
"""
function getxy(a::Array, N::Int64=80)
    x = cosine_spacing(N)
    y = zeros(N)
    if real(a)!=a
        x = complex(x)
        y = complex(y)
    end
    for i=1:N
        for j=1:6
            y[i] += a[j]*x[i]^(j-0.5)
        end #for coeffs
    end #for x stations

    return y
end

"""
    solveas(p::Array, up::Bool, directTE::Bool)

Calculate the PARSEC coefficients using modified parameters (see parsec() docstring) for either the top or bottom curve.

### parameters:
[1]: rLE
[2]: X
[3]: Z
[4]: Zxx
[5]: θTE
[6]: ZTE
"""
function solveas(p::Array, up::Bool)
    #---define b in Ax=b
        if up==true
            b = [sqrt(2*p[1]); p[3]; 0.0; p[4]; tan(p[5]); p[6]]
        else
            b = [-sqrt(2*p[1]); p[3]; 0.0; p[4]; tan(p[5]); p[6]]
        end
    #---define A in Ax=b
    #rename for convenience
    x = p[2]
    #calculate exponenets and coefficients
    n2 = zeros(6)
    n3 = zeros(6)
    n4 = zeros(6)
    c3 = zeros(6)
    c4 = zeros(6)
    for i=1:6
        n2[i] = i-0.5
        n3[i] = i-1.5
        n4[i] = i-2.5
        c3[i] = i-0.5
        c4[i] = (i-0.5)*(i-1.5)
    end
    c6 = n2'

    #assemble matrix rows
    r1 = [ 1.0 0.0 0.0 0.0 0.0 0.0 ]
    r2 = [ x^n2[1] x^n2[2] x^n2[3] x^n2[4] x^n2[5] x^n2[6] ]
    r3 = [ c3[1]*x^n3[1] c3[2]*x^n3[2] c3[3]*x^n3[3] c3[4]*x^n3[4] c3[5]*x^n3[5] c3[6]*x^n3[6] ]
    r4 = [ c4[1]*x^n4[1] c4[2]*x^n4[2] c4[3]*x^n4[3] c4[4]*x^n4[4] c4[5]*x^n4[5] c4[6]*x^n4[6] ]
    r5 = c6
    r6 = ones(1,6)

    #assemble matrix
    A = [ r1; r2; r3; r4; r5; r6 ]

    #solve for x in Ax=b (A\b)
    return A\b
end

"""
    parsec(p::Array; N::Int=80)

Calculate the x,y airfoil coordinates for both top and bottom surfaces using modified PARSEC Parameterization method.

Use parsecstd() for standard parsec implementation.  This modified version employs direct values for trailing edge position and angles for each surface.

### parameters:
[1] - rLE
[2] - Xup
[3] - Xlo
[4] - Zup
[5] - Zlo
[6] - Zxxup
[7] - Zxxlo
[8] - θTEup
[9] - θTElo
[10] - ZTEup (optional)
[11] - ZTElo (optional)
[12] - AoA (permitted, unused)

### keyword arguments:
N: Number of x stations along chord
"""
function parsec(p::Array, N::Int64=80)
    if length(p) == 9
        p = [p; 0.0; 0.0] #ZTEup = ZTElo = 0
    elseif length(p) == 12
        p = p[1:end-1]
    elseif length(p) != 11
        error("Incorrect number of parameters, must have one of the following: \n9: Sharp TE \n11: Full Parameter Set \n12: Full Set + AoA")
    end
    #--- Get x-values ---#
    x = cosine_spacing(N)

    #--- Upper Curve ---#
    au = solveas([ p[1]; p[2]; p[4]; p[6]; p[8]; p[10]], true)
    yu = getxy(au, N)

    #--- Lower Curve ---#
    al = solveas([ p[1]; p[3]; p[5]; p[7]; p[9]; p[11] ], false)
    yl = getxy(al, N)

    return x, yu, yl
end

"""
    xytoparsec(x::Array{Float64,1}, yup::Array{Float64,1}, ylo::Array{Float64,1})

Calculate approximate modified (see parsec() docstring) PARSEC parameters based on input x,y coordinates.

**The following parameters are not accurate:
- Xlo
- Zlo
- Zxxlo
- θTEup
- θTElo
(as of 28 Jan 2019)**
"""
function xytoparsec(x::Array{Float64,1}, yup::Array{Float64,1}, ylo::Array{Float64,1})
    #spline curves
    spup = Dierckx.Spline1D(x,yup)
    splo = Dierckx.Spline1D(x,ylo)
    h = 1e-6

    # dzle = spup(x[1]+h/2) - splo(x[1]-h/2)
    # dzzle = (spup(x[1]+h) - (splo(x[1])+spup(x[1])) + splo(x[1]-h))/(h^2)
    # rLE = abs(((1+dzle^2)^3/2)/dzzle)
    #-- LE radius --#
    #get points at LE
    if isapprox(x[1],0.0,atol=1e-15)
        p = [x[2]; ylo[2]]
        q = [x[1]; ylo[1]]
        r = [x[2]; yup[2]]
    else
        p = [x[1]; ylo[1]]
        q = [0.0; 0.0]
        r = [x[1]; yup[1]]
    end
    #define function to find circle centerpoint
    function findC!(C,x)
        C[1] = (p[1]-x[1])^2 + (p[2]-x[2])^2 - (q[1]-x[1])^2 - (q[2]-x[2])^2
        C[2] = (q[1]-x[1])^2 + (q[2]-x[2])^2 - (r[1]-x[1])^2 - (r[2]-x[2])^2
    end
    #solve for center point
    center = NLsolve.nlsolve(findC!,[1.0; 1.0])
    #calculate radius
    rLE = sqrt((q[1]-center.zero[1])^2 + (q[2]-center.zero[2])^2)

    #-- TE Position --#
    if isapprox(yup[end] ,0.0, atol=1e-15)
        ZTEup = 0.0
    else
        ZTEup = yup[end]
    end

    if isapprox(ylo[end] ,0.0, atol=1e-15)
        ZTElo = 0.0
    else
        ZTElo = ylo[end]
    end

    #-- TE Angles --#
    #lower surface TE angle
    thetaTElo = atan(ylo[end-1]-ylo[end],x[end-1]-x[end])

    #upper surface TE angle
    thetaTEup = atan(yup[end-1]-yup[end],x[end-1]-x[end])

    #-- Max Thicknesses and Related Curvatures --#

    #maximum y coordinate on top surface
    Zup, idxup = findmax(yup)
    #x coordinate for maximum y coordinate
    Xup = x[idxup]
    #average panel length on either side of max y coordinate
    # hup = mean([x[idxup]-x[idxup-1]; x[idxup+1]-x[idxup]])

    #curvature at minimum y coordinate (2nd order central diff)
    Zxxup = (spup(x[idxup]+h) - 2*spup(x[idxup]) + spup(x[idxup]-h))/(h^2)

    #minimum y coordinate on bottom surface
    Zlo, idxlo = findmin(ylo)
    #x coordinate for minimum y coordinate
    Xlo = x[idxlo]
    # #average panel length on either side of max y coordinate
    # hlo = mean([x[idxlo-1]-x[idxlo]; x[idxlo]-x[idxlo+1]])

    #curvature at minimum y coordinate (2nd order central diff)
    Zxxlo = (splo(x[idxlo]+h) - 2*splo(x[idxlo]) + splo(x[idxlo]-h))/(h^2)

    return[rLE; Xup; Xlo; Zup; Zlo; Zxxup; Zxxlo; thetaTEup; thetaTElo; ZTEup; ZTElo]
end

"""
    xytoparsecLsqFit(x,zu,zl)

Uses LsqFit to go from xy to modified parsec without calculating each parameter.
Depending on the initial guess and subsequent iterations, it may throw an error
involving complex numbers. If so, alter the initial guess.
"""
function xytoparsecLsqFit(x,zu,zl)

    function model(xx,p)
        x = xx[1:div(end,2)] #splits x into 2 vectors for my function
        xi,zu,zl = parsec(p,length(x))
        return vcat(zu,zl)
    end

    #initial guess
    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

    fit = LsqFit.curve_fit(model,vcat(x,x),vcat(zu,zl),guess)
    p = fit.param'

    return p
end

"""
    solveasstd(p::Array, up::Bool, directTE::Bool)

Calculate the PARSEC coefficients using standard parameters (see parsec() docstring) for either the top or bottom curve.

### parameters:
[1]: rLE
[2]: X
[3]: Z
[4]: Zxx
[5]: αTE
[6]: βTE
[7]: ΔZTE
[8]: ZTE
"""
function solveasstd(p::Array, up::Bool)
    #---define b in Ax=b
        if up==true
            b = [sqrt(2*p[1]); p[8]+p[7]/2.0; p[3]; tan(p[5]-p[6]); 0.0; p[4]]
        else
            b = [-sqrt(2*p[1]); p[8]+p[7]/2.0; p[3]; tan(p[5]-p[6]); 0.0; p[4]]
        end
    #---define A in Ax=b
    #rename for convenience
    x = p[2]
    #calculate exponenets and coefficients
    n2 = zeros(6)
    n3 = zeros(6)
    n4 = zeros(6)
    c3 = zeros(6)
    c4 = zeros(6)
    for i=1:6
        n2[i] = i-0.5
        n3[i] = i-1.5
        n4[i] = i-2.5
        c3[i] = i-0.5
        c4[i] = (i-0.5)*(i-1.5)
    end
    c6 = n2'

    #assemble matrix rows
    r1 = [ 1.0 0.0 0.0 0.0 0.0 0.0 ]
    r2 = ones(1,6)
    r3 = [ x^n2[1] x^n2[2] x^n2[3] x^n2[4] x^n2[5] x^n2[6] ]
    r4 = c6
    r5 = [ c3[1]*x^n3[1] c3[2]*x^n3[2] c3[3]*x^n3[3] c3[4]*x^n3[4] c3[5]*x^n3[5] c3[6]*x^n3[6] ]
    r6 = [ c4[1]*x^n4[1] c4[2]*x^n4[2] c4[3]*x^n4[3] c4[4]*x^n4[4] c4[5]*x^n4[5] c4[6]*x^n4[6] ]

    #assemble matrix
    A = [ r1; r2; r3; r4; r5; r6 ]

    #solve for x in Ax=b (A\b)
    return A\b
end

"""
    parsecstd(p::Array, N::Int64=80)

Calculate the x,y airfoil coordinates for both top and bottom surfaces using standard PARSEC Parameterization method.

Use parsec() for modified parsec implementation.

### parameters:
[1] - rLE
[2] - Xup
[3] - Xlo
[4] - Zup
[5] - Zlo
[6] - Zxxup
[7] - Zxxlo
[8] - αTE
[9] - βTE
[10] - ΔZTE (optional)
[11] - ZTE (optional)
[12] - AoA (permitted, unused)

### keyword arguments:
N: Number of x stations along chord
"""
function parsecstd(p::Array, N::Int64=80)
    if length(p) == 9
        p = [p; 0.0; 0.0] #ZTE = ΔZTE = 0
    elseif length(p) == 12
        p = p[1:end-1]
    elseif length(p) != 11
        error("Incorrect number of parameters, must have one of the following: \n9: Sharp TE \n11: Full Parameter Set \n12: Full Set + AoA")
    end
    #--- Get x-values ---#
    x = cosine_spacing(N)

    #--- Upper Curve ---#
    au = solveasstd([ p[1]; p[2]; p[4]; p[6]; p[8]; p[9]; p[10]; p[11]], true)
    yu = getxy(au, N)

    #--- Lower Curve ---#
    al = solveasstd([ p[1]; p[3]; p[5]; p[7]; p[8]; p[9]; p[10]; p[11]], false)
    yl = getxy(al, N)

    return x, yu, yl
end

#"""
#    xytoparsecLsqFitstd(x,zu,zl)

#Uses LsqFit to go from xy to traditional parsec instead of calculating each parameter.
#Depending on the initial guess and subsequent iterations, it may throw an error
#involving imaginary numbers. If so, alter the initial guess.
#"""
#function xytoparsecLsqFitstd(x,zu,zl)

#    function model(xx,p)
#        x = xx[1:div(end,2)] #splits x into 2 vectors for my function
#        xi,zu,zl = parsecstd(p,length(x))
#        return vcat(zu,zl)
#    end

#    #initial guess
#    guess = [0.01, 0.5, 0.5, 0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.0]

#    fit = LsqFit.curve_fit(model,vcat(x,x),vcat(zu,zl),guess)
#    return fit

#end

#"""
#Not Complete
#"""
#function xytoparsecstd()

#    #-- TE Angles --#
#    #y coord of mean camber line 1 station forward from TE
#    yc1 = (yup[end-1]-ylo[end-1])/2
#    #vector of camberline
#    vcam = [x[end]-x[end-1]; yup[end]-yc1]
#    #vector of chordline
#    vcord = [-1.0; 0.0]
#    #angle of mean camber line at TE
#    alphaTE = -acos(dot(vcam,vcord)/(norm(vcam)*norm(vcord)))

#    #panel vector of upper surface adjacent to TE
#    vup = [x[end-1]-x[end]; yup[end-1]-yup[end]]
#    #panel vector of lower surface adjacent to TE
#    vlo = [x[end-1]-x[end]; ylo[end-1]-ylo[end]]
#    #Angle between upper and lower surfaces at TE
#    betaTE = acos(dot(vup,vlo)/(norm(vup)*norm(vlo)))

#end
