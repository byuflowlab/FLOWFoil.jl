#---------------------------------#
#          Reformat Inputs        #
#---------------------------------#

"""
    reformat_inputs(x, y, flow_angles)
    reformat_inputs(coordinates, flow_angles)

Reformats inputs to be the expected format.

Specifically, reverses coordinates that appear to be backwards, and puts coordinates in the vector of matrix format used throughout the various methods.

# Arguments
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `flow_angles::Vector{Float}` : Vector of angles of attack in degrees (may be a single float as well)
OR
- `x::Vector{Float}` : Vector of x-coordinates of airfoil geometry
- `y::Vector{Float}` : Vector of y-coordinates of airfoil geometry
- `flow_angles::Vector{Float}` : Vector of angles of attack in degrees (may be a single float as well)

# Returns
- `coordinates::Vector{Float}` : reformatted coordinates
- `nbody::Int` : number of bodies
- `aoa::Vector{Float}` : reformatted angles of attack
"""
function reformat_inputs(x, y, flow_angles)
    # combine coordinates
    return reformat_inputs([x y], flow_angles)
end

function reformat_inputs(coordinates, flow_angles)
    # if coordinates is a tuple of coordinates, then splat it into a vector
    return reformat_inputs([coordinates...], flow_angles)
end

function reformat_inputs(coordinates::AbstractMatrix, flow_angles)
    sc = size(coordinates)
    @assert length(sc) == 2 "Coordinates must include both x and y coordinates"
    @assert 2 ∈ sc "Coordinates must only include x and y values"

    if sc[1] == 2
        # shape into expected dimensions
        return reformat_inputs([coordinates[1, :] coordinates[2, :]], flow_angles)
    else
        # pass through if already in expected dimensions
        return reformat_inputs([coordinates], flow_angles)
    end
end

function reformat_inputs(coordinates::AbstractArray, flow_angles)

    # - Get Number of Bodies - #
    if length(size(coordinates)) == 1
        # if coordinates is a vector of vectors
        # should be if user inputs it as such or if user inputs a tuple
        nbody = length(coordinates)
    else
        # if coordinates is a single matrix
        # should be if user inputs a single airfoil or if separate x and y coordinates are combined
        nbody = 1
    end

    # - Make Angles of Attack a Vector (if not one already) - #
    if length(flow_angles) == 1
        aoa = [flow_angles[1]]
    else
        aoa = flow_angles
    end

    # - Check that ordering is correct - #
    for ni in 1:nbody
        if abs(coordinates[ni][1, 1] - coordinates[ni][end, 1]) > 1e-1
            if coordinates[ni][1, 1] > coordinates[ni][end, 1]
                @info "It appears body $(ni) is a body of revolution and the coordinates are given in reverse order.  Reversing coordinates to be from front to back. If this is not a body of revolution, then your trailing edge gap is likely too large."
                reverse!(coordinates[ni]; dims=1)
            end
        else
            nc = ceil(Int, length(coordinates[ni][:, 1]) / 2)
            if sum(coordinates[ni][1:nc, 2]) / nc > sum(coordinates[ni][nc:end, 2]) / nc
                @info "It appears that body $(ni) is an airfoil with coordinates given in reverse order. Reversing the coordinates to be clockwise. If the coordinates are in the correct order, then your airfoil is likely upside down."
                reverse!(coordinates[ni]; dims=1)
            end
        end
    end

    return coordinates, nbody, aoa#, reynolds_numbers, mach_numbers
end

#---------------------------------#
#          Panel Geometry         #
#---------------------------------#

"""
    generate_panel_geometry(p, coordinates)

Generate panel geometries for a given set of coordinates.

# Arguments:
- `method::Method` : method for solving
- `coordinates::Vector{Matrix{Float}}` : Vector of [x y] matrices of airfoil coordinates (may be a single matrix as well)

# Returns:
- `panel_geometry::Vector{NTuple}` : Vector of named tuples (one for each body in the system)
"""
function generate_panel_geometry(method::Method, coordinates) end

#---------------------------------#
#         System Geometry         #
#---------------------------------#

"""
    generate_system_geometry(method::Method, panel_geometry; gap_tolerance=1e-10) end

Generate relative geometry between panels used in assembling the system matrices.

# Arguments:
- `method::Method` : Problem type object for dispatch
- `panel_geometry::Vector{Panel}` : Array of panel object for airfoil system. (can also be a single panel object if only one body is being modeled)

# Returns:
- `system_goemetry::NTuple` : Named tuple including various influence geometries for the system
- `TE_geometry::NTuple` : Named tuple specifically for trailing edge gap panel_geometry if present.
"""
function generate_system_geometry(method::Method, panel_geometry; gap_tolerance=1e-10) end

#---------------------------------#
#   Generate (Non)Linear System   #
#---------------------------------#

"""
    generate_system_matrices(method::Method, mesh, TEmesh)

Assemble various matrices used in system solves.

# Arguments:
- `method::Method` : method object for dispatch
- `panel_geometry::Vector{NTuple}` : Vector of named tuples (one for each body in the system)
- `system_geometry::NTuple` : geometry for airfoil system to analyze.
- `TE_geometry::NTuple` : Trailing edge gap panel influence geometry.

# Returns:
- `system::NTuple` : Named tuple containing relevant influence and right hand side matrices, etc.
"""
function generate_system_matrices(
    method::Method, panel_geometry, system_geometry, TE_geometry
) end

#---------------------------------#
#              SOLVE              #
#---------------------------------#

"""
    solve(method::Method, system_matrices)

Solve system.

# Arguments:
- `method::Method` : method object for dispatch
- `system_matrices::system_matrices` : system to solve

# Returns:
 - `strengths:Array{Float}` : values for solved panel strengths
"""
function solve(method::Method, system_matrices) end

#---------------------------------#
#          Post Processing        #
#---------------------------------#

"""
    post_process(method::Method, panel_geometry, system_geometry, strengths)

Post-process solution and produce a Polar object.

# Arguments:
- `method::Method` : Problem type for dispatch
- `panel_geometry::Vector{Panel}` : vector of Panel objects
- `system_geometry::NTuple` : geometry for airfoil system to analyze.
- `strengths:Array{Float}` : values for solved panel strengths

# Returns:
- `outputs::Outputs` : object of type outputs
"""
function post_process(method::Method, panel_geometry, system_geometry, strengths) end
