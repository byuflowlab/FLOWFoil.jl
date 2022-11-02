#=
Singularity Distributions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
10/22 - Add axisymmetric ring vortex and related functions
=#

"""
    get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\overline{\\Psi}^\\gamma\$

**Arguments:**
 - `theta1::Float` : angle between panel and vector from node1 to evaluation point
 - `theta2::Float` : angle between panel and vector from node2 to evaluation point
 - `ln1::Float` : value of ln(rmag1), which may be that or 0.0, depening on evaluation point location
 - `ln2::Float` : value of ln(rmag2), which may be that or 0.0, depening on evaluation point location
 - `dmag::Float` : panel length
 - `h::Float` : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - `a::Float` : length of base of right triangle with height, h, and hypontenuse, r1.

"""
function get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)
    return 1.0 / (2.0 * pi) * (h * (theta2 - theta1) - dmag + a * ln1 - (a - dmag) * ln2)
end

"""
    get_psitildegamma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\widetilde{\\Psi}^\\gamma\$

**Arguments:**
 - `psibargamma::Float` : value of \$\\overline{\\Psi}^\\gamma\$
 - `r1mag::Float` : distance from node1 to evaluation point
 - `r2mag::Float` : distance from node2 to evaluation point
 - `theta1::Float` : angle between panel and vector from node1 to evaluation point
 - `theta2::Float` : angle between panel and vector from node2 to evaluation point
 - `ln1::Float` : value of ln(rmag1), which may be that or 0.0, depening on evaluation point location
 - `ln2::Float` : value of ln(rmag2), which may be that or 0.0, depening on evaluation point location
 - `dmag::Float` : panel length
 - `h::Float` : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - `a::Float` : length of base of right triangle with height, h, and hypontenuse, r1.

"""
function get_psitildegamma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)
    ptg =
        a * psibargamma +
        1.0 / (4.0 * pi) * (r2mag^2 * ln2 - r1mag^2 * ln1 - r2mag^2 / 2.0 + r1mag^2 / 2.0)
    return (dmag == 0.0) ? 0.0 : (ptg / dmag)
end

"""
    get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\overline{\\Psi}^\\sigma\$

**Arguments:**
 - `theta1::Float` : Angle between panel and evaluation point, centered at node1.
 - `theta2::Float` : Angle between panel and evaluation point, centered at node2.
 - `ln1::Float` : Natural log of distance from node1 to evaluation point.
 - `ln2::Float` : Natural log of distance from node2 to evaluation point.
 - `h::Float` : Distance from panel to evaluation in panel normal direction.
 - `a::Float` : Distance from node1 to evaluation in panel tangent direction.
"""
function get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)
    return 1 / (2 * pi) * (a * (theta1 - theta2) + dmag * theta2 + h * ln1 - h * ln2)
end

"""
    get_psitildesigma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\widetilde{\\Psi}^\\sigma\$

**Arguments:**
 - `psibargamma::Float` : value of \$\\overline{\\Psi}^\\sigma\$
 - `r1mag::Float` : distance from node1 to evaluation point
 - `r2mag::Float` : distance from node2 to evaluation point
 - `theta1::Float` : angle between panel and vector from node1 to evaluation point
 - `theta2::Float` : angle between panel and vector from node2 to evaluation point
 - `dmag::Float` : panel length
 - `h::Float` : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - `a::Float` : length of base of right triangle with height, h, and hypontenuse, r1.
"""
function get_psitildesigma(psibarsigma, r1mag, r2mag, theta1, theta2, dmag, h, a)
    return a / dmag * phibarsigma +
           1.0 / (4 * pi * dmag) * (rmag2^2 * theta2 - rmag1^2 * theta1 - h * dmag)
end

"""
    get_vortex_influence(node1, node2, point)

Calculate vortex influence coefficients on the evaluation point from the panel between node1 and node2.

**Arguments:**
 - `node1::Array{Float}(2)` : [x y] location of node1
 - `node2::Array{Float}(2)` : [x y] location of node2
 - `point::Array{Float}(2)` : [x y] location of evaluation point

"""
function get_vortex_influence(node1, node2, point)

    # Use inputs to get raw distances
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Calculate a, h, and natural logs based on position of point
    theta1, theta2, ln1, ln2, h, a = get_orientation(node1, node2, point)

    # get psibargamma value
    psibargamma = get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)

    # get psitildegamma value
    psitildegamma = get_psitildegamma(
        psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a
    )
    # put psi`s together
    return (psibargamma - psitildegamma), psitildegamma
end

"""
    get_source_influence(node1, node2, point)

Calculate source influence coefficients on the evaluation point from the panel between node1 and node2.

**Arguments:**
 - `node1::Array{Float}(2)` : [x y] location of node1
 - `node2::Array{Float}(2)` : [x y] location of node2
 - `point::Array{Float}(2)` : [x y] location of evaluation point
"""
function get_source_influence(node1, node2, point)

    # Use inputs to get raw distances
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Calculate a, h, and natural logs based on position of point
    theta1, theta2, ln1, ln2, h, a = get_orientation(node1, node2, point)

    #get psibarsigma value
    psibarsigma = get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)

    # shift source in order to get a better behaved branch cut orientation
    if (theta1 + theta2) > pi
        psibarsigma -= 0.25 * dmag
    else
        psibarsigma += 0.75 * dmag
    end

    return psibarsigma
end

####################################
##### ----- AXISYMMETRIC ----- #####
####################################

"""
    get_ring_vortex_influence(paneli, panelj)

Cacluate the influence of a ring vortex at panel j onto panel i.

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

**Returns:**
- `aij::Float` : Influence of vortex ring strength at panel j onto panel i.
"""
function get_ring_vortex_influence(paneli, panelj)

    #get geometry of panels and influence
    x, r, rj, dmagj, m, nhati = get_ring_geometry(paneli, panelj)

    #calculate unit velocities
    u = get_u_ring(x, r, rj, dmagj, m)
    v = get_v_ring(x, r, rj, m)

    #return appropriate strength
    # if asin(sqrt(m)) != pi / 2
    if m != 1.0

        #panels are different
        return (-u * cos(paneli.beta) + v * sin(paneli.beta)) * dmagj
    else
        #same panel -> self induction equation

        #NOTE: this is not eqn 4.22 in Lewis.  Their code uses this expression which seems to avoid singularities better.  Not sure how they changed the second term (from dj/4piR to -R) though; perhaps the R in the text != the curvature in the code (radiusofcurvature vs curvature).

        # constant used in multiple places to clean things up
        cons = 4.0 * pi * rj / dmagj

        # return self inducement coefficient
        return -0.5 - panelj.radiusofcurvature -
               (log(2.0 * cons) - 0.25) / cons * cos(panelj.beta)
    end
end

"""
    get_u_ring(x, r, rj, m)

Calculate x-component of velocity influence of vortex ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `uij::Float` : x-component of velocity induced by panel j onto panel i
"""
function get_u_ring(x, r, rj, dj, m; probe=false)

    #get the first denominator
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * (r - 1)
    den2 = x^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #phi = asin(sqrt(m))
    #if probe && phi > 89.5 * pi / 180.0
    #    println("K, E: ", K, ", ", E)
    #end
    ##return velocity
    #if sqrt((x^2 + (r - 1.0)^2)) < 0.01 && probe
    #    println("u sing: ", (r - 1.0) / (2.0 * pi * sqrt(x^2 + (r - 1.0)^2)))
    #    return (r - 1.0) / (2.0 * pi * sqrt(x^2 + (r - 1.0)^2))
    #else
    return 1.0 / den1 * (K - (1.0 + num2 / den2) * E)
    # end
end

"""
    get_v_ring(x, r, rj, m)

Calculate r-component of velocity influence of vortex ring.

**Arguments:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

**Returns:**
- `vij::Float` : r-component of velocity induced by panel j onto panel i
"""
function get_v_ring(x, r, rj, m; probe=false)

    #get numerator and denominator of first fraction
    num1 = x / r
    den1 = 2.0 * pi * rj * sqrt(x^2 + (r + 1.0)^2)

    num2 = 2 * r
    den2 = x^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #return velocity
    # if sqrt((x^2 + (r - 1.0)^2)) < 0.01 && probe
    #     return -x / (2.0 * pi * sqrt(x^2 + (r - 1.0)^2))
    # else
    return num1 / den1 * (K - (1.0 + num2 / den2) * E)
    # end
end

"""
    get_elliptics(m)

Calculate value of elliptic functions for the given geometry parameter.

**Arguments:**
- `m::Float` : Elliptic Function parameter

**Returns:**
- `K::Float` : K(m), value of elliptic function of the first kind at m.
- `E::Float` : E(m), value of eeliptic function of the second kind at m.
"""
function get_elliptics(m)

    #calculate phi
    # phi = asin(sqrt(m))

    #if phi > 89.5 * pi / 180.0
    #    #if singular, use asymptotic expressions
    #    K = log(4.0 / cos(phi))
    #    E = 1.0 + 0.5 * (K - 1.0 / 1.2) * cos(phi)^2
    #    return K, E
    #else
    #looks like special functions uses some sort of asymptotic or equivalent expressions already.
    if m > 1 || isnan(m)
        #m cannot be greater than 1 for elliptic functions, and cannot mathematically be either, but numerically might be infinitesimally larger.
        m = 1.0
    end
    return SpecialFunctions.ellipk(m), SpecialFunctions.ellipe(m)
    # end
end
