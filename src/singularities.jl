#=
Singularity Distributions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
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

