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

"""
"""
function generate_inviscid_system(problem::Problem)
    return generate_inviscid_system(problem.mesh, problem.solver)
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
- `mesh::Array{PlanarMesh}` : PlanarMesh for airfoil to analyze.
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
