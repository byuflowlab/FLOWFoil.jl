#=

Convenience functions wrapping problem, system, solution, and post processing steps into single functions for user convenience

Authors: Judd Mehr,

=#
"""
    outputs(cl, cd, cm, cp, vs)

# Returns:
- `outputs::Struct` : Struct with outputs.  Nominally contains
  - `cl`: lift coefficient of each body
  - `cd`: total drag coefficient of each body
  - `cm`: moment coefficient of each body
  - `cp`: pressure coefficient for each panel of each body
  - `vs`: surface velocities on each body
"""
struct outputs{TF}
    cl::Matrix{TF}
    cd::Matrix{TF}
    cm::Matrix{TF}
    cp::Vector{Matrix{TF}}
    vs::Vector{Matrix{TF}}
end

"""
    analyze(coordinates, flow_angles=0.0, reynolds=1e6, machs=0.0; method::Method=Mfoil())
    analyze(x, y, flow_angles=0.0, reynolds=1e6, machs=0.0; method::Method=Mfoil())

Convenience function for setting up, solving, and post-processing airfoils and airfoil systems.

# Required Arguments:
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
OR
- `x::Vector{Float}` : Vector of x-coordinates of airfoil geometry
- `y::Vector{Float}` : Vector of y-coordinates of airfoil geometry

Note that inputting separate vectors for airfoil coordinates is only available for analysis of single airfoils/bodies.  Multi-airfoil/body systems require the use of a tuple of matrices for coordinate inputs.

# Optional Arguments:
- `flow_angles::Vector{Float}=0.0` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}=1e-6` : Vector of reynolds numbers (may be a single float as well)
- `machs::Vector{Float}=0.0` : Vector of machs numbers (may be a single float as well)

Note that Reynolds and Mach numbers are only used for viscous methods, and Flow Angles are unused in the axisymmetric methods.

# Keyword Arguments:
- `method::Method=Mfoil()` : desired method for solving

# Returns:
- `outputs::NTuple` : named tuple with outputs.  Nominally contains
  - `cl`: lift coefficient of each body
  - `cd`: total drag coefficient of each body
  - `cdp`: profile drag coefficient of each body
  - `cm`: moment coefficient of each body
  - `tangential_velocities`: surface velocities on each body
  - `surface_pressures`: surface pressures on each body
  - `convergenced`: convergence flag
  - `auxiliary outputs`: a named tuple that contains additional outputs applicable to the method used.
"""
function analyze(
    x::AbstractVector,
    y::AbstractVector,
    flow_angles=[0.0],
    reynolds=[1e6],
    machs=[0.0];
    method::Method=Mfoil()
)
    return analyze([x y], flow_angles, reynolds, machs; method=method)
end

function analyze(
    coordinates, flow_angles=[0.0], reynolds=[1e6], machs=[0.0]; method::Method=Mfoil()
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
    return post_process(method, panel_geometry, system_geometry, strengths, flow_angles)
end
