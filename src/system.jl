#=

Inviscid System Functions

Authors: Judd Mehr,

=#

"""
    InviscidSystem

**Fields:**
 - `A::Array{Float,2}` : Coefficient Matrix on Left Hand Side.
 - `b::Array{Float,2}` : Boundary Condition Coefficient Vector on Right Hand Side.
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
struct InviscidSystem{TA,TB,TI}
    A::TA
    b::TB
    Ns::TI
end

######################################################################
#                                                                    #
#                         SYSTEM GENERATION                          #
#                                                                    #
######################################################################
"""
"""
function generate_inviscid_system(problem::Problem)
    return generate_inviscid_system(problem.mesh, problem.type)
end

"""
**Arguments:**
- `mesh::Array{PlanarMesh}` : PlanarMesh for airfoil to analyze.
"""
function generate_inviscid_system(mesh, ::Planar)
    # Get coeffiecient matrix (A, left hand side)
    A, Ns = assemble_ring_vortex_matrix(mesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_ring_boundary_conditions(mesh)

    return InviscidSystem(A, b, Ns)
end

"""
**Arguments:**
- `mesh::Array{AxisymMesh}` : AxisymMesh for airfoil to analyze.
"""
function generate_inviscid_system(mesh, ::Axisymmetric)
    # Get coeffiecient matrix (A, left hand side)
    A, Ns = assemble_vortex_matrix(mesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_boundary_conditions(mesh)

    return InviscidSystem(A, b, Ns)
end

"""
**Arguments:**
- `mesh::Array{PlanarMesh}` : PlanarMesh for airfoil to analyze.
"""
function generate_inviscid_system(mesh, ::Periodic) end

######################################################################
#                                                                    #
#                         UTILITY FUNCTIONS                          #
#                                                                    #
######################################################################

"""
    size_system(mesh)

Count size of inviscid system matrix.

**Arguments:**
 - `mesh::Array{Mesh}` : The system for which to calculate the linear system size.
"""
function size_system(meshes, ::Planar)

    # initialize
    # number of bodies for convenience
    numbodies = length(meshes)

    # initialize total system size
    N = 0

    # initialize system size contributions from each mesh
    Ns = [1 for i in 1:numbodies]

    # Count number of airfoil nodes in each mesh.
    for i in 1:numbodies
        Ns[i] = length(meshes[i].nodes)
        N += Ns[i]
    end

    return N, Ns
end

function size_system(meshes, ::Axisymmetric)

    # initialize
    # number of bodies for convenience
    numbodies = length(meshes)

    # initialize total system size
    N = 0

    # initialize system size contributions from each mesh
    Ns = [1 for i in 1:numbodies]

    # Count number of airfoil nodes in each mesh.
    for i in 1:numbodies
        Ns[i] = length(meshes[i].panels)
        N += Ns[i]
    end

    return N, Ns
end

"""
    get_offset(Ns)

Get the offset values for the mesh system to be used in the system matrix assembly.

**Arguments:**
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
function get_offset(Ns)
    return [0; cumsum(Ns[1:(end - 1)])]
end
