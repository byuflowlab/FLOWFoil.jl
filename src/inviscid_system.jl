#=
Inviscid System Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
10/22 - Add axisymmetric solution capabilities
=#

"""
    get_inviscid_system(meshes)

Calculate, then gather the vortex and boundary condition matricies into an InviscidSystem object.

**Arguments:**
- `meshes::Array{PlanarMesh}` : PlanarMesh for airfoil to analyze.
"""
function get_inviscid_system(meshes; axisymmetric=false)
    if !axisymmetric
        # Get coeffiecient matrix (A, left hand side)
        Aij, Ns = assemble_vortex_matrix(meshes)

        # Get boundary conditions (RHS)
        b = assemble_boundary_conditions(meshes)
    else
        # get coefficient matrix
        Aij, Ns = assemble_ring_vortex_matrix(meshes)

        # get boundary conditions
        b = assemble_ring_boundary_conditions(meshes)
    end

    return InviscidSystem(Aij, b, Ns)
end

####################################
##### ----- AXISYMMETRIC ----- #####
####################################
