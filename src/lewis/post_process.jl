# NOTE: flow_angles is unused, just here to make convenience functions work everywhere
function post_process(method::Lewis, panel_geometry, system_geometry, strengths, flow_angles; npanels=80)
    return post_process(
        method, [panel_geometry], system_geometry, strengths, flow_angles; npanels=npanels
    )
end

function post_process(
    method::Lewis, panel_geometry::AbstractVector, system_geometry, strengths, flow_angles; npanels=80
)

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    nbodies = system_geometry.nbodies

    # - Initialize Outputs - #
    TF = eltype(system_geometry.k2)

    vs = [zeros(idx[m][end]-idx[m][1]+1, 1) for m in 1:nbodies]
    cp = [zeros(idx[m][end]-idx[m][1]+1, 1) for m in 1:nbodies]
    cl = zeros(1, nbodies)
    cd = zeros(1, nbodies)
    cm = zeros(1, nbodies)

    for m in 1:nbodies
        # - Extract surface velocity - #
        #= Note that we here assume that we are using the subtractive method for the kutta conditions, requiring us to recover the strengths value for the last panel (trailing edge upper side).  We also assume here that the indexing starts at the lower side trailing edge and proceeds clockwise back to the upper side trailing edge.  Otherwise, not only will the solver not have worked, there will also be an indexing error here
        =#
        vs[m][:, 1] = [
            strengths[idx[m][1:(end - 1)]]
            -strengths[idx[m][1]]
        ]

        # - Calculate surface pressure - #
        cp[m][:, 1] = 1.0 .- (vs[m][:, 1]) .^ 2
    end

    if nbodies == 1
        return InviscidOutputs(vs[1], cp[1], cl, cd, cm)
    else
        return InviscidOutputs(vs, cp, cl, cd, cm)
    end
end

"""
    calculate_duct_thrust(inviscid_solution; Vinf=1.0, rho=1.225)

Calculate the thrust of the duct.
TODO: need to test!

# Arguments:

# Keyword Arguments:
- `Vinf::Float` : freestream velocity
- `rho::Float` : air density

# Returns:
- `thrust::Float` : duct thrust (negative indicates drag)
"""
function calculate_duct_thrust(
    method::Lewis, outputs, panel_geometry, system_geometry; Vinf=1.0, rho=1.225
)

    # - Rename for Convenience - #
    idx = system_geometry.node_indices

    # Calculate dynamic pressure
    q = 0.5 * rho * Vinf^2

    # Initialize output
    TF = eltype(system_geometry.m)
    fx = zeros(TF, system_geometry.nbodies)

    # Loop through Bodies
    for m in 1:(system_geometry.nbodies)
        if !method.body_of_revolution[m]

            #dimensionalize pressure
            P = outputs.surface_pressure[idx[m]] .* q

            # add panel pressure in x-direction to total sectional force
            fx += sum(
                P .* panel_geometry[m].panel_length[:] .*
                panel_geometry[m].panel_normal[:, 1] .*
                panel_geometry[m].center_point[:, 2],
            )
        end
    end

    #return total duct thrust for whole annulus: -fx*2pi
    return fx * 2.0 * pi
end
