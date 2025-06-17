"""
    InviscidOutputs

Note: not all methods return values for all outputs.  Methods will return zeros in such cases.

# Fields
- `vs`: surface velocities normalized by freestream velocity on each body, nominally a matrix, but becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cp`: pressure coefficient for each panel of each body, becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cl`: lift coefficient of each body, nominally a vector but becomes a matrix in the multi-body case with dimensions angle x body
- `cd`: Inviscid drag coefficient of each body (simply integral of pressure coefficient in x direction), but becomes a matrix in the multi-body case with dimensions angle x body
- `cm`: moment coefficient of each body, but becomes a matrix in the multi-body case with dimensions angle x body
"""
struct InviscidOutputs{TM,TV}
    vs::TM
    cp::TM
    cl::TV
    cd::TV
    cm::TV
end

"""
    analyze(coordinates, flow_angles=0.0; method::Method=Mfoil())
    analyze(x, y, flow_angles=0.0; method::Method=Mfoil())

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
- `outputs::OutputType` : outputs object (note that only inviscid methods are currently implemented)
"""
function analyze(x, y, flow_angles; method::Method=Mfoil())
    return analyze([x y], flow_angles; method=method)
end

function analyze(coordinates, flow_angles; method::Method=Mfoil())

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

        # Post Process Solution
        return post_process(method, panel_geometry, system_geometry, strengths, flow_angles)
    end
end
