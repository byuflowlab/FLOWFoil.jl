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

"""
    assemble_ring_vortex_matrix(meshes)

Assemble vortex coefficient matrix with full N+nk x N+nk system, including Kutta condition(s), where nk is the number of required kutta conditions in the system, and N is the total number of nodes between all the meshes.

**Arguments:**
 - `meshes::Array{AxiSymMesh}` : Mesh System for which to solve.

**Returns:**
 - `amat::Matrix{Float,2}` : aij (LHS) coefficient matrix including kutta condition augmentation if required.
"""
function assemble_ring_vortex_matrix(meshes)

    # size sysetm
    N, Ns = FLOWFoil.size_system(meshes; axisymmetric=true)
    n = length(Ns)
    nk = countkutta(meshes)
    offset = FLOWFoil.get_offset(Ns)

    # initialize coefficient matrix
    amat = [0.0 for i in 1:(N + nk), j in 1:(N + nk)]

    # Loop through system
    for x in 1:n
        for y in 1:n
            if x == y && !meshes[x].bodyofrevolution
                #call the kutta version of coefiicient matrix to do back substitution step

                amat[(1 + offset[x]):(Ns[x] + offset[x]), (1 + offset[y]):(Ns[y] + offset[y])] .= FLOWFoil.assemble_ring_vortex_coefficients(
                    meshes[x], meshes[y]; backsub=true
                )

                # put in the kutta condition for each airfoil (end rows of the system matrix)
                amat[N + x, 1 + offset[y]] = 1.0
                amat[N + x, Ns[y] + offset[y]] = 1.0

                #put unit bound vortex value in each row
                amat[(1 + offset[x]):(Ns[x] + offset[x]), N + y] .= 1.0

            else

                # get influence coefficients for each portion of the sysetm (mesh y acts on mesh x)
                # put things in the correct place in the system matrix
                amat[(1 + offset[x]):(Ns[x] + offset[x]), (1 + offset[y]):(Ns[y] + offset[y])] .= FLOWFoil.assemble_ring_vortex_coefficients(
                    meshes[x], meshes[y]
                )
            end
        end
    end

    return amat, Ns
end

"""
    assemble_ring_vortex_coefficients(meshi, meshj; kutta)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxM portion of the system influence coefficient matrix associated with the M-1 panels of meshj acting on the N nodes of meshi. It does not include the kutta condition or the influence of the constant stream function on the airfoil nodes.

**Arguments:**
 - `meshi::AxiSymMesh` : mesh being influenced.
 - `meshj::AxiSymMesh` : mesh doing the influencing.

**Keyword Arguments:**
- `backsub:Bool` : flag as to whether a backsubstitution is to be applied (happens for annular airfoils when meshi==meshj)

**Returns:**
- `amat::Matrix{Float,2}` : NxN matrix of aij coefficients (does not include kutta condition augmentation).
"""
function assemble_ring_vortex_coefficients(meshi, meshj; backsub=false)

    # get nodes for convenience
    panelsi = meshi.panels
    panelsj = meshj.panels

    # get system size
    N = length(panelsi)
    M = length(panelsj)

    #initialize NxN coefficient matrix
    amat = [0.0 for i in 1:N, j in 1:M]

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:M

            # obtain influence coefficient for ith evaluation point and jth panel
            amat[i, j] = FLOWFoil.get_ring_vortex_influence(panelsi[i], panelsj[j])
        end
    end

    if backsub

        #apply back substitution to matrix
        for i in 1:N
            sum = 0.0
            jidx = N + 1 - i
            for j in 1:M
                if j != jidx
                    sum += amat[j, i] * meshj.panels[j].length
                end
            end
            dmagj = meshj.panels[jidx].length
            amat[jidx, i] = -sum / dmagj
        end

        # TODO: doesn't seem to change anything...
        # # Bound Vortex Correction
        # for i in 1:N
        #     for j in 1:M
        #         amat[i, j] += meshj.panels[j].length
        #     end
        # end

        return amat
    else
        return amat
    end
end

"""
    assemble_ring_boundary_conditions(meshes)

Assemble boundary condition (RHS) vector.

**Arguments:**
 - `meshes::Array{AxiSymMesh}` : mesh system for which to solve

**Returns**
 - `bc::Array{Float}` : Boundary condition array.
"""
function assemble_ring_boundary_conditions(meshes)

    # size the system
    N, Ns = size_system(meshes; axisymmetric=true)
    nk = countkutta(meshes)
    offset = get_offset(Ns)

    # initialize boundary condition array
    bc = [0.0 for i in 1:(N + nk)]

    # Loop through system
    for m in 1:length(Ns)

        # get node locations for convenience
        panels = meshes[m].panels
        N = length(panels)

        # generate boundary condition array
        if meshes[m].bodyofrevolution
            bc[(1 + offset[m]):(Ns[m] + offset[m]), 1] = [-cos(panels[i].beta) for i in 1:N]
        else
            bc[(1 + offset[m]):(Ns[m] + offset[m]), 1] = [
                1.0 - cos(panels[i].beta) for i in 1:N
            ]
        end
    end

    return bc
end

"""
    countkutta(meshes)

Count the number of meshes for which to apply the kutta condition.

**Arguments:**
- `meshes::Array{AxiSymMesh}` : mesh system over which to count.

**Returns:**
- `nk::Int` : number of meshes requiring kutta condition treatment.
"""
function countkutta(meshes)

    #count how many meshes have kutta condtion
    nk = 0
    for i in 1:length(meshes)
        nk += meshes[i].bodyofrevolution ? 0 : 1
    end

    return nk
end
