#=

Convenience functions wrapping problem, system, solution, and post processing steps into single functions for user convenience

Authors: Judd Mehr, Ayden Bennett

=#

"""
    InviscidOutputs

Note: not all methods return values for all outputs.  Methods will return zeros in such cases.

# Fields
- `vs`: surface velocities on each body, nominally a matrix, but becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cp`: pressure coefficient for each panel of each body, becomes a vector of matrices in the multi-body case with dimensions [body][panel,angle]
- `cl`: lift coefficient of each body, nominally a vector but becomes a matrix in the multi-body case with dimensions angle x body
- `cd`: total drag coefficient of each body, but becomes a matrix in the multi-body case with dimensions angle x body
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
    analyze(coordinates, flow_angles=0.0, reynolds=1e6, mach=0.0; method::Method=Mfoil())
    analyze(x, y, flow_angles=0.0, reynolds=1e6, mach=0.0; method::Method=Mfoil())

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Required Arguments
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
OR
- `x::Vector{Float}` : Vector of x-coordinates of airfoil geometry
- `y::Vector{Float}` : Vector of y-coordinates of airfoil geometry

Note that inputting separate vectors for airfoil coordinates is only available for analysis of single airfoils/bodies.  Multi-airfoil/body systems require the use of a tuple of matrices for coordinate inputs.

# Optional Arguments
- `flow_angles::Vector{Float}=0.0` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}=1e-6` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}=0.0` : Vector of mach numbers (may be a single float as well)

Note that Reynolds and Mach numbers are only used for viscous methods, and Flow Angles are unused in the axisymmetric methods.

# Keyword Arguments
- `method::Method=Mfoil()` : desired method for solving

# Returns
- `outputs::InviscidOutputs` : outputs object (note that only inviscid methods are currently implemented)
"""
function analyze(x, y, flow_angles; method::Method=Mfoil())
    return analyze([x y], flow_angles; method=method)
end

function analyze(coordinates, flow_angles; method::Method=Mfoil())

    # Reformat inputs as needed
    coordinates, nbodies, flow_angles = reformat_inputs(coordinates, flow_angles)

    if typeof(method) <: NeuralFoil
        return analyze_nf(coordinates[1], flow_angles; method=method)
    else
    if typeof(method) <: Martensen
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
end