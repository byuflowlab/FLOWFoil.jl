function post_process(
    ap::AxisymmetricProblem, problem, panels::TP, mesh, solution; npanels=80, debug=false
) where {TP<:Panel}
    return post_process(
        ap::AxisymmetricProblem, problem, [panels], mesh, solution; npanels=80, debug=false
    )
end

function post_process(
    ap::AxisymmetricProblem, problem, panels, mesh, solution; npanels=80, debug=false
)

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    nbodies = mesh.nbodies

    # - Initialize Outputs - #
    TF = eltype(mesh.m)
    v_surf = zeros(TF, nbodies, 2 * npanels - 1)
    p_surf = zeros(TF, nbodies, 2 * npanels - 1)
    xsmooth = zeros(TF, nbodies, 2 * npanels - 1)

    for m in 1:nbodies
        # - Extract surface velocity - #
        #= Note that we here assume that we are using the subtractive method for the kutta conditions, requiring us to recover the solution value for the last panel (trailing edge upper side).  We also assume here that the indexing starts at the lower side trailing edge and proceeds clockwise back to the upper side trailing edge.  Otherwise, not only will the solver not have worked, there will also be an indexing error here
        =#
        vti = [solution.x[idx[m][1:(end - 1)]]; -solution.x[idx[m][1]]]

        # - Calculate surface pressure - #
        cpi = 1.0 .- (vti) .^ 2

        ### --- Smooth Distributions --- ###
        #smooth_distributions functions are found in utils.jl
        v_surf[m, :], xsmooth[m, :] = smooth_distributions(
            Constant(),
            panels[m].panel_center,
            vti,
            npanels;
            body_of_revolution=ap.body_of_revolution[m],
        )
        p_surf[m, :], _ = smooth_distributions(
            Constant(),
            panels[m].panel_center,
            cpi,
            npanels;
            body_of_revolution=ap.body_of_revolution[m],
        )
    end

    return AxisymmetricPolar(v_surf, p_surf, xsmooth)
end

"""
    calculate_duct_thrust(inviscid_solution; Vinf=1.0, rho=1.225)

Calculate the thrust of the duct.
TODO: need to test!

**Arguments:**
- `inviscid_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `Vinf::Float` : freestream velocity

**Keyword Arguments:**
- `rho::Float` : air density

**Returns:**
- `thrust::Float` : duct thrust (negative indicates drag)
"""
function calculate_duct_thrust(polar::AxisymmetricPolar, panels, mesh; Vinf=1.0, rho=1.225)

    # - Rename for Convenience - #
    idx = mesh.node_indices

    # Calculate dynamic pressure
    q = 0.5 * rho * Vinf^2

    # Initialize output
    TF = eltype(mesh.m)
    fx = zeros(TF, mesh.nbodies)

    # Loop through Bodies
    for m in 1:(mesh.nbodies)
        if !problem.body_of_revolution[m]

            #dimensionalize pressure
            P = polar.surface_pressure[idx[m]] .* q

            # add panel pressure in x-direction to total sectional force
            fx += sum(
                P .* panels[m].panel_length[:] .* panels[m].panel_normal[:, 1] .*
                panels[m].center_point[:, 2],
            )
        end
    end

    #return total duct thrust for whole annulus: -fx*2pi
    return fx * 2.0 * pi
end
