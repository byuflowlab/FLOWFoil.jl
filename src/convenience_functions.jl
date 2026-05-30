"""
    InviscidOutputs

Note: not all methods return values for all outputs.  Methods will return zeros in such cases.

# Fields
- `vs`: surface velocities normalized by freestream velocity on each body, nominally a matrix, but becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cp`: pressure coefficient for each panel of each body, becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cl`: lift coefficient of each body, nominally a vector but becomes a matrix in the multi-body case with dimensions angle x body
- `cd`: Inviscid drag coefficient of each body (simply integral of pressure coefficient in x direction), but becomes a matrix in the multi-body case with dimensions angle x body
- `cm`: moment coefficient of each body, but becomes a matrix in the multi-body case with dimensions angle x body
- `gamma_ref`: (Xfoil only) raw node vortex strengths from the linear solve, organized as a (nnodes × 2) matrix whose columns hold the strengths for unit cos α and unit sin α respectively. `nothing` for methods that don't expose it.
- `A`: (Xfoil only) the inviscid linear-system influence matrix used to obtain `gamma_ref`. `nothing` for methods that don't expose it.
"""
struct InviscidOutputs{TM,TV,TG,TA}
    vs::TM
    cp::TM
    cl::TV
    cd::TV
    cm::TV
    gamma_ref::TG
    A::TA
end

InviscidOutputs(vs, cp, cl, cd, cm) = InviscidOutputs(vs, cp, cl, cd, cm, nothing, nothing)

"""
    analyze(coordinates, flow_angles=0.0; method::Method=Xfoil())
    analyze(x, y, flow_angles=0.0; method::Method=Xfoil())

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Arguments
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `flow_angles::Vector{Float}` : Vector of angles of attack in degrees (may be a single float as well)
OR
- `x::Vector{Float}` : Vector of x-coordinates of airfoil geometry
- `y::Vector{Float}` : Vector of y-coordinates of airfoil geometry
- `flow_angles::Vector{Float}` : Vector of angles of attack in degrees (may be a single float as well)

Note that inputting separate vectors for airfoil coordinates is only available for analysis of single airfoils/bodies.  Multi-airfoil/body systems require the use of a tuple of matrices for coordinate inputs.

# Keyword Arguments
- `method::MethodType` : desired method for solving

# Returns
- For inviscid methods (default `Xfoil()`, `Lewis()`, `Martensen()`, etc.):
  an `InviscidOutputs` struct with fields `vs, cp, cl, cd, cm, gamma_ref, A`.
- For viscous Xfoil (`Xfoil(viscous=true, Re=…, …)`): a flat NamedTuple
  with both the inviscid intermediates (`gamma_ref`, `A`) and the per-α
  viscous results (`cl, cd, cm, transition_upper, transition_lower,
  converged, per_alpha`). The viscous path requires a single body and a
  blunt-trailing-edge airfoil; multi-α is supported by looping internally.
"""
function analyze(x, y, flow_angles; method::Method=Xfoil())
    return analyze([x y], flow_angles; method=method)
end

function analyze(coordinates, flow_angles; method::Method=Xfoil())

    # Reformat inputs as needed
    coordinates, nbodies, flow_angles = reformat_inputs(coordinates, flow_angles)

    if typeof(method) <: NeuralFoil
        return analyze_nf(coordinates[1], flow_angles; method=method)
    elseif typeof(method) <: LegacyXfoil
        return analyze_lxf(coordinates[1], flow_angles; method=method)
    else
        # Generate Panel Geometry
        panel_geometry = generate_panel_geometry(method, coordinates)

        # Generate Influence Mesh
        system_geometry = generate_system_geometry(method, panel_geometry)

        # Assemble Linear System
        system_matrices = generate_system_matrices(method, panel_geometry, system_geometry)

        # Solve System
        strengths = solve(method, system_matrices)

        # Dispatch viscous coupled solve (Xfoil method only) or inviscid post-processing
        if typeof(method) <: Xfoil && method.viscous
            return analyze_viscous(
                method,
                panel_geometry,
                system_geometry,
                system_matrices,
                strengths,
                flow_angles,
            )
        else
            return post_process(
                method,
                panel_geometry,
                system_geometry,
                system_matrices,
                strengths,
                flow_angles,
            )
        end
    end
end
