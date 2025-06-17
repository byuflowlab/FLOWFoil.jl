"""
    calculate_ring_vortex_influence(paneli, panelj, system_geometry, i, j)

Cacluate the influence of a ring vortex at panel j onto panel i.

# Arguments
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `system_geometry::FLOWFoil.AxisymmetricMesh` : relative geometry object.
- `i::Int` : index for ith panel
- `j::Int` : index for jth panel

# Returns
- `aij::Float` : Influence of vortex ring strength at panel j onto panel i.
"""
function calculate_ring_vortex_influence(paneli, panelj, system_geometry, i, j)
    m2p = system_geometry.mesh2panel

    #calculate unit velocities
    u = get_u_ring_vortex(
        system_geometry.y[i, j],
        system_geometry.r[i, j],
        panelj.panel_center[m2p[j], 2],
        panelj.panel_length[m2p[j]],
        system_geometry.k2[i, j],
    )

    v = get_v_ring_vortex(
        system_geometry.y[i, j],
        system_geometry.r[i, j],
        panelj.panel_center[m2p[j], 2],
        system_geometry.k2[i, j],
    )

    #return appropriate strength
    # if asin(sqrt(m)) != pi / 2
    if system_geometry.k2[i, j] != 1.0

        #panels are different
        return (u * cos(paneli.panel_angle[m2p[i]]) + v * sin(paneli.panel_angle[m2p[i]])) *
               panelj.panel_length[m2p[j]]
    else
        #same panel -> self induction equation

        #NOTE: this is not eqn 4.22 in Lewis.  Their code uses this expression which seems to avoid singularities better.  Not sure how they changed the second term (from dj/4piR to -R) though; perhaps the R in the text != the curvature in the code (radiusofcurvature vs curvature).

        # constant used in multiple places to clean things up
        cons = 4.0 * pi * panelj.panel_center[m2p[j], 2] / panelj.panel_length[m2p[j]]

        # return self inducement coefficient
        return -0.5 - panelj.panel_curvature[m2p[j]] -
               (log(2.0 * cons) - 0.25) / cons * cos(panelj.panel_angle[m2p[j]])
    end
end

"""
    get_u_ring_vortex(y, r, rj, m)

Calculate y-component of velocity influence of vortex ring.

# Arguments
- `y::Float` : ratio of difference of ith and jth panel y-locations and jth panel r-location ( (zi-zj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

# Returns
- `uij::Float` : y-component of velocity induced by panel j onto panel i
"""
function get_u_ring_vortex(y, r, rj, dj, m; probe=false)

    #get the first denominator
    den1 = 2.0 * pi * rj * sqrt(y^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * (r - 1)
    den2 = y^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    return -1.0 / den1 * (K - (1.0 + num2 / den2) * E)
end

"""
    get_v_ring_vortex(y, r, rj, m)

Calculate r-component of velocity influence of vortex ring.

# Arguments
- `y::Float` : ratio of difference of ith and jth panel y-locations and jth panel r-location ( (zi-zj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

# Returns
- `vij::Float` : r-component of velocity induced by panel j onto panel i
"""
function get_v_ring_vortex(y, r, rj, m; probe=false)

    #get numerator and denominator of first fraction
    num1 = y / r
    den1 = 2.0 * pi * rj * sqrt(y^2 + (r + 1.0)^2)

    num2 = 2 * r
    den2 = y^2 + (r - 1)^2

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    return num1 / den1 * (K - (1.0 + num2 / den2) * E)
end

"""
    get_elliptics(m)

Calculate value of elliptic functions for the given geometry parameter.

# Arguments:
- `m::Float` : Elliptic Function parameter

# Returns
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

"""
    get_u_ring_source(y, r, rj, m)

Calculate y-component of velocity influence of source ring.

# Arguments
- `y::Float` : ratio of difference of ith and jth panel y-locations and jth panel r-location ( (zi-zj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

# Returns
- `uij::Float` : y-component of velocity induced by panel j onto panel i
"""
function get_u_ring_source(y, r, rj, dj, m)
    #TODO: dj unused, remove from inputs and throughout

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #get the first denominator
    den1 = 2.0 * pi * rj * sqrt(y^2 + (r + 1.0)^2)

    #get numerator and denominator of second fraction
    num2 = 2 * y * E
    den2 = y^2 + (r - 1)^2

    return 1.0 / den1 * (num2 / den2)
end

"""
    get_v_ring_source(y, r, rj, m)

Calculate r-component of velocity influence of source ring.

# Arguments
- `y::Float` : ratio of difference of ith and jth panel y-locations and jth panel r-location ( (zi-zj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `m::Float` : Elliptic Function parameter

# Returns
- `vij::Float` : r-component of velocity induced by panel j onto panel i
"""
function get_v_ring_source(y, r, rj, m)

    #get values for elliptic integrals
    K, E = get_elliptics(m)

    #get numerator and denominator of first fraction
    den1 = 2.0 * pi * rj * sqrt(y^2 + (r + 1.0)^2)

    num2 = 2 * r * (r - 1.0)
    den2 = y^2 + (r - 1)^2

    return 1.0 / den1 * (K - (1.0 - num2 / den2) * E)
end
