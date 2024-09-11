#=

Convenience functions wrapping problem, system, solution, and post processing steps into single functions for user convenience

Authors: Judd Mehr,

=#

"""
    analyze(coordinates, flow_angle, reynolds, mach, method)

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Arguments:
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)

# Optional Arguments:
- `flow_angle::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)

# Keyword Arguments:
- `method::Method=Xfoil()` : desired method for solving
- `gap_tolerance::Float` : gap_tolerance for determining trailing edge gap

# Returns:
- `outputs::Outputs` : object of type Outputs
"""
function analyze(
        x, y, flow_angle=[0.0], reynolds=[1e6], mach=[0.0]; method::Mfoil(), gap_tolerance=1e-10
)
    return analyze(
        [x y], flow_angle, reynolds, mach; method=method, gap_tolerance=gap_tolerance
    )
end

function analyze(
    coordinates,
    flow_angle=[0.0],
    reynolds=[1e6],
    mach=[0.0];
    method::Mfoil(),
    gap_tolerance=1e-10,
)

    # Reformat inputs as needed
    problem = reformat_inputs(coordinates, flow_angle, reynolds, mach)

    # Generate Panel Geometry
    panel_geometry = generate_panel_geometry(method, coordinates)

    # Generate Influence Mesh
    system_geometry, TE_geometry = generate_system_geometry(
        method, panel_geometry; gap_tolerance=gap_tolerance
    )

    # Assemble Linear System
    system_matrices = generate_system_matrices(
        method, panel_geometry, system_geometry, TE_geometry
    )

    # Solve System
    strengths = solve(method, system_matrices)

    # Post Process Solution
    return post_process(method, panel_geometry, system_geometry, strengths)
end
