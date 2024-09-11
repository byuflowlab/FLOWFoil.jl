#---------------------------------#
#             Methods             #
#---------------------------------#

abstract type Method end

#---------------------------------#
#          Reformat Inputs        #
#---------------------------------#

"""
    reformat_inputs(method::method, coordinates, flow_angle, reynolds, mach)

Defines a Problem object to be used for setup and post-processing.

**Arguments:**
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `flow_angle::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)
- `method::Method` : Type of problem to solve (planar, axisymmetric, or periodic).

**Returns:**
- `problem::Problem` : Problem object
"""
function reformat_inputs(method::Method, coordinates, flow_angle, reynolds, mach)

    # - Get Number of Bodies - #
    #= NOTE:
        The size function in the first if will use the overloaded version that returns the length of the tuple if coordinates is of type Tuple.  See overloaded function at the end of FLOWFoil.jl
    =#
    if length(size(coordinates)) == 1
        nbody = length(coordinates)
    else
        nbody = size(coordinates)[2]
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

    # - Must be Viscous if Reynolds Numbers are Greater Than Zero - #
    if all(x -> x != nothing, reynolds_numbers) && all(x -> x > 0.0, reynolds_numbers)
        viscous = true
    elseif all(x -> x == nothing, reynolds_numbers) || all(x -> x <= 0.0, reynolds_numbers)
        viscous = false
    else
        @error "Cannot mix viscid and inviscid analyses. For viscous analsis sets, all reynolds numbers must be greater than 0.0. For inviscid analysis sets, it is best to set reynolds to -1.0."
    end

    # - RETURN - #
    return coordinates, nbody, aoa, reynolds_numbers, mach_numbers
end
#---------------------------------#
#          Panel Geometry         #
#---------------------------------#

"""
    generate_panels(p, coordinates)

Generate panel object for a give set of coordinates.

**Arguments:**
- `p::Method` : problem type object
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)

**Returns:**
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
"""
function generate_panels(::Method, coordinates) end

#---------------------------------#
#         System Geometry         #
#---------------------------------#
"""
**Arguments:**
- `method::Method` : Problem type object for dispatch
- `panels::Vector{Panel}` : Array of panel object for airfoil system. (can also be a single panel object if only one body is being modeled)

**Returns:**
- `mesh::Mesh` : Mesh object including various influence geometries for the system
- `TEmesh::Mesh` : Mesh object specifically for trailing edge gap panels if present.
"""
function generate_system_geometry(method::Method, panels; gap_tolerance=1e-10) end

#---------------------------------#
#   Generate (Non)Linear System   #
#---------------------------------#

"""
    generate_system_matrices(method::Method, mesh, TEmesh)

**Arguments:**
- If PlanarProblem:
  - `method::Method` : method object for dispatch
  - `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
  - `mesh::Mesh` : Mesh for airfoil system to analyze.
  - `TEMesh::Mesh` : Trailing edge gap panel influence mesh.

- If AxisymmetricProblem:
  - `method::Method` : method object for dispatch
  - `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
  - `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
  - `mesh::Mesh` : Mesh for airfoil system to analyze.

**Returns:**
` inviscid_system::InviscidSystem` : Inviscid System object containing influence and boundary condition matrices for the system.
"""
function generate_system_matrices(method::Method, panels, mesh, TEmesh) end

#---------------------------------#
#              SOLVE              #
#---------------------------------#

"""
    solve(problem)

Solve problem defined by the input Problem object and return the solution in a Solution object.

**Arguments:**
- `problem::Problem` : Problem to solve

**Returns:**
 - `solution::{InviscidSolution or ViscousSolution}` : returns solution of type matching viscous flag in problem.
"""
function solve(system)
    solution = solve_inviscid(system)

    return solution
end

#---------------------------------#
#          Post Processing        #
#---------------------------------#

abstract type Outputs end
abstract type AuxOutputs end

"""
    post_process(::Method, problem, panels, mesh, solution; debug=false)

Post-process solution and produce a Polar object.

**Arguments:**
- `method::Method` : Problem type for dispatch
- `problem::Problem` : Problem object
- `panels::Vector{Panel}` : vector of Panel objects
- `mesh::Mesh` : Mesh object
- `solution::Solution` : Solution object

**Keyword Arguments:**
- `npanels::Int` : number of panels to use on top and bottom surface for surface distribution smoothing
- `debug::Bool` : Flag for output format (see below)

**Returns:**
- If debug == false: `polar::Polar` : a Polar object
- If debug == true: all the fields of a Polar object (in order), but as a tuple rather than a struct, such that the raw, unsmoothed, velocity and pressure distributions can be output.
"""
function post_process(::Method, problem, panels, mesh, solution; npanels=80, debug=false) end
