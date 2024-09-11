#---------------------------------#
#             Methods             #
#---------------------------------#

abstract type Method end

#---------------------------------#
#          Reformat Inputs        #
#---------------------------------#

"""
    reformat_inputs(x, y, flow_angle, reynolds, mach)
    reformat_inputs(coordinates, flow_angle, reynolds, mach)

Reformats inputs to be the expected format.

# Arguments:
- `coordinates::Tuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a vecotr of matrices as well)
- `flow_angle::Vector{Float}` : Vector of angles of attack (may be a single float as well)

# Optional Arguments:
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)

# Keyword Arguments:
- `method::Method=Xfoil()` : desired method for solving
- `gap_tolerance::Float` : gap_tolerance for determining trailing edge gap

# Returns:
- `coordinates::Vector{Float}` : reformatted coordinates
- `nbody::Int` : number of bodies
- `aoa::Vector{Float}` : reformatted angles of attack
- `reynolds_numbers::Vector{Float}` : reformatted Reynolds numbers
- `mach_numbers::Vector{Float}` : reformatted Mach numbers
"""
function reformat_inputs(x, y, flow_angle, reynolds, mach)
    # combine coordinates
    return reformat_inputs([x y], flow_angle, reynolds, mach)
end

function reformat_inputs(coordinates, flow_angle, reynolds, mach)
    # if coordinates is a tuple of coordinates, then splat it into a vector
    return reformat_inputs([coordinates...], flow_angle, reynolds, mach)
end

function reformat_inputs(coordinates::AbstractMatrix, flow_angle, reynolds, mach)
    sc = size(coordinates)
    @assert length(sc) == 2 "Coordinates must include both x and y coordinates"
    @assert 2 ∈ sc "Coordinates must only include x and y values"

    if sc[1] == 2
        # shape into expected dimensions
        return reformat_inputs(
            [coordinates[1, :] coordiantes[2, :]], flow_angle, reynolds, mach
        )
    else
        # pass through if already in expected dimensions
        return reformat_inputs(coordinates, flow_angle, reynolds, mach)
    end
end

function reformat_inputs(coordinates::AbstractArray, flow_angle, reynolds, mach)

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
    if length(flow_angle) == 1
        aoa = [flow_angle[1]]
    else
        aoa = flow_angle
    end

    # - Make Reynolds Number a Vector (if not one already) - #
    if reynolds == nothing
        reynolds_numbers = [reynolds[1]]
    elseif length(reynolds) == 1
        reynolds_numbers = [reynolds[1]]
    else
        reynolds_numbers = reynolds
    end

    # - Make Mach Number a Vector (if not one already) - #
    if mach == nothing
        mach_numbers = [mach[1]]
    elseif length(mach) == 1
        mach_numbers = [mach[1]]
    else
        mach_numbers = mach
    end

    return coordinates, nbody, aoa, reynolds_numbers, mach_numbers
end

#---------------------------------#
#          Panel Geometry         #
#---------------------------------#

"""
    generate_panel_geometry(p, coordinates)

Generate panel geometries for a give set of coordinates.

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
function solve(method::Method, system) end

#---------------------------------#
#          Post Processing        #
#---------------------------------#

abstract type Outputs end
abstract type AuxOutputs end

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
