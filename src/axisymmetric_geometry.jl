
"""
    AxiSymMesh{TP,TB}

Axisymmetric Mesh Object

**Fields:**
- `panels::FLOWFoil.AxiSymPanel` : panel objects describing surface geometry.
- `bodyofrevolution::Bool` : Flag as to whether or not the mesh represents a body of revolution.
"""
struct AxiSymMesh{TP,TB}
    panels::TP
    bodyofrevolution::TB
end

"""
    AxiSymPanel{TF,TA}

Panel object for axisymmetric meshes.

**Fields:**
- `controlpoint::Array{Float}` : [x;r] coordinates of panel midpoint.
- `length::Float` : length of panel
- `normal::Array{Float}` : unit normal vector of panel (TODO: remove if unused)
- `beta::Float` : angle panel makes with positive x-axis (radians)
- `radiusofcurvature::Float` : the radius of curvature of the geometry at the panel control point. TODO: make sure this is actually correct with current implementation.
"""
struct AxiSymPanel{TF,TA}
    controlpoint::TA
    length::TF
    normal::TA
    beta::TF
    radiusofcurvature::TF
end
