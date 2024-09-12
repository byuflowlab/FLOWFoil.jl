#=

Convenience functions wrapping problem, system, solution, and post processing steps into single functions for user convenience

Authors: Judd Mehr,

=#

"""
    analyze(coordinates, flow_angles, reynolds, machs, method)

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Arguments:
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)

# Optional Arguments:
- `flow_angles::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `machs::Vector{Float}` : Vector of machs numbers (may be a single float as well)

# Keyword Arguments:
- `method::Method=Xfoil()` : desired method for solving

# Returns:
- `outputs::Outputs` : object of type Outputs
"""
function analyze(
    x::AbstractVector,
    y::AbstractVector,
    flow_angles=[0.0],
    reynolds=[1e6],
    machs=[0.0];
    method=Mfoil(),
    gap_tolerance=1e-10,
)
    return analyze(
        [x y], flow_angles, reynolds, machs; method=method, gap_tolerance=gap_tolerance
    )
end

function analyze(
    coordinates,
    flow_angles=[0.0],
    reynolds=[1e6],
    machs=[0.0];
    method=Mfoil(),
    gap_tolerance=1e-10,
)

    # Reformat inputs as needed
    coordinates, nbodies, flow_angles, reynolds, machs = reformat_inputs(
        coordinates, flow_angles, reynolds, machs
    )

    # Generate Panel Geometry
    panel_geometry = generate_panel_geometry(method, coordinates)

    # Generate Influence Mesh
    system_geometry = generate_system_geometry(method, panel_geometry)

    # Assemble Linear System
    system_matrices = generate_system_matrices(method, panel_geometry, system_geometry)

    # Solve System
    strengths = solve(method, system_matrices)

    # Post Process Solution
    return post_process(method, panel_geometry, system_geometry, strengths)
end
