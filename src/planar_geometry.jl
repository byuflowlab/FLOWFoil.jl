
"""
    PlanarMesh{TF,TB,TN}

Mesh for single body.

**Fields:**
 - `nodes::Array{Array{Float,2}}` : [x y] node (panel edge) locations for airfoil
 - `chord::Float` : airfoil chord length
 - `blunt_te::Bool` : boolean for whether or not the trailing edge is blunt or not.
 - `trailing_edge_gap::Float` : trailing edge gap distance
 - `tdp::Float` : dot product of unit vectors of trailing edge bisection and gap vectors
 - `txp::Float` : pseudo-cross product of unit vectors of trailing edge bisection and gap vectors
**Assuptions:**
 - x and y coordinates start at the bottom trailing edge and proceed clockwise.

"""
struct PlanarMesh{TF,TB,TN<:Vector{Matrix{TF}}}
    nodes::TN
    chord::TF
    blunt_te::TB
    trailing_edge_gap::TF
    tdp::TF
    txp::TF
end

"""
    PlanarMeshSystem{TM,TF,TL}

System of meshes to solve.

**Fields:**
 - `meshes::Array{Mesh}` : Array of mesh objects.
 - `scales::Vector{Float}` : Airfoil scaling factors.
 - `angles::Vector{Float}` : Airfoil angles of attack.
 - `locations::Array{Array{TF}}` : Array of leading edge locations.

"""
struct PlanarMeshSystem{TM,TF,TL<:Vector{Matrix{TF}}}
    meshes::TM
    scales::TF
    angles::TF
    locations::TL
end
