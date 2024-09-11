"""
    calculate_periodic_vortex_influence(paneli, panelj)

Cacluate the influence of a periodic vortex at panel j onto panel i.
TODO: BUGGY, NOT WORKING, NEED TO FIX

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

**Returns:**
- `aij::Float` : Influence of vortex strength at panel j onto panel i.
"""
function calculate_periodic_vortex_influence(::Constant, paneli, panelj, mesh, i, j)

    # - Self Induction Term - #
    if isapprox([mesh.x[i, j]; mesh.y[i, j]], [0.0; 0.0])
        return -0.5 - paneli.delta_angle[i] / (4.0 * pi)

    else

        # - Standard Periodic Coefficient - #
        s = mesh.stagger * pi / 180.0
        t = mesh.pitch

        a = (mesh.x[i, j] * cos(s) - mesh.y[i, j] * sin(s)) * 2 * pi / t
        b = (mesh.x[i, j] * sin(s) + mesh.y[i, j] * cos(s)) * 2 * pi / t
        e = exp(a)
        sinha = 0.5 * (e - 1.0 / e)
        cosha = 0.5 * (e + 1.0 / e)
        k = 0.5 / t / (cosha - cos(b))

        dmagi = paneli.panel_length[i]
        dmagj = panelj.panel_length[j]

        betaj = panelj.panel_angle[j]

        return sinha * sin(s + betaj) - sin(b) * cos(s + betaj) * k * dmagj
    end
end
