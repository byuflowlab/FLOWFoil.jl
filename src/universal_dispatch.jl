#---------------------------------#
#             Methods             #
#---------------------------------#

abstract type Method end

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
