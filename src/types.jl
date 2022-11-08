#=
Compose Type Definitions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
10/22 - Axisymmetric types added
=#

"""
    Problem{TM,TF,TB}

Problem definition (geometry, operating point(s), and method selection) and output behavior.

**Fields:**
 - `meshes::Array{PlanarMesh}` : Array of mesh objects
 - `angleofattack::Float` : angle of attack to analyze.
 - `reynolds::Float` : Reynolds number to analyze.
 - `mach::Float` : Mach number to analyze.
 - `viscous::Bool` : Flag to solve viscous or inviscid only
 - `verbose::Bool` : Flag to print out verbose statements
"""
struct Problem{TM,TF,TB}
    meshes::TM
    angleofattack::TF
    reynolds::TF
    mach::TF
    viscous::TB
    axisymmetric::TB
    verbose::TB
end

"""
    Problem(meshes, angleofattack=0.0, reynolds=0.0, mach=0.0; viscous=true, verbose=false)

Constructor for Problem Objects.

**Arguments:**
 - `meshes::Array{PlanarMesh or AxiSymMesh}` : Array of mesh objects
 - `angleofattack::Float` : Angle of Attack (currently unused)
 - `reynolds::Float` : Reynolds Number (currently unused)
 - `mach::Float` : Mach Number (currently unused)

**Keyword Arguments:**
 - `viscous::Bool` : Flag to solve viscous or inviscid only
 - `axisymmetric::Bool` : Flag for axisymmetric solver.
 - `verbose::Bool` : Flag to print out verbose statements
"""
function Problem(
    meshes::TM,
    angleofattack::TF1=0.0,
    reynolds::TF2=0.0,
    mach::TF3=0.0;
    viscous::TB=true,
    axisymmetric::TB=false,
    verbose::TB=false,
) where {TM, TF1, TF2, TF3, TB}

    TF = promote_type(TF1, TF2, TF3)

    return Problem{TM, TF, TB}(meshes, angleofattack, reynolds, mach, viscous, axisymmetric, verbose)
end

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

# """
#     WakeMesh{TF}

# **Fields:**
#  - `wake_nodes::Array{Float,2}` : x,y wake panel node locations.
#  - `wake_midpoints::Array{Float,2}` : x,y wake panel center point locations.
# """
# struct WakeMesh{TF}
#     wake_nodes::Array{TF}
#     wake_midpoints::Array{TF}
# end

"""
    InviscidSystem{TF}

**Fields:**
 - `vcoeffmat::Array{Float,2}` : Vortex Coefficient Matrix used in solution.
 - `bccoeffvec::Array{Float,2}` : Boundary Coefficient Vector used in solution.
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
struct InviscidSystem{TA,TB,TI}
    vcoeffmat::TA
    bccoeffvec::TB
    Ns::TI
end

# """
# **Fields:**
# """
# struct ViscousSystem{TPa,TPr,TF}
#     parameters::TPa
#     properties::TPr
#     vcoeffmat::Array{TF}
#     scoeffmat::Array{TF}
#     bccoeffvec::Array{TF}
# end

"""
    InviscidSolution{TM,TF,TD}

**Fields:**
 - `mesh::PlanarMesh` : PlanarMesh object describing airfoil nodes etc.
 - `panelgammas::Array{Float,2}` : \$\\gamma_0\$ and \$\\gamma_{90}\$ values at each airfoil node.
 - `bodystrength::Array{Float}` : if 2D system, bodystrength = \$\\Psi_0\$ (constant stream function) 0 and 90 values.  If axisymmetric system, bodystrength = bound vortex strength of body.
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
 - `system::InviscidSystem` : system object.
"""
struct InviscidSolution{TM,TF,TI,TD}
    meshes::TM
    panelgammas::TF
    bodystrength::TF
    Ns::TI
    system::TD
end

"""
    Polar{TF}

**Fields:**
 - `lift::Float` : Lift Coefficient.
 - `drag::Float` : Total Drag Coefficient.
 - `pdrag::Float` : Pressure Drag Coefficient.
 - `idrag::Float` : Induced Drag Coefficient.
 - `moment::Float` : Moment Coefficient.
 - `surfacevelocity::Vector{Float}` : surface velocity distribution
 - `surfacepressure::Vector{Float}` : surface pressure distribution
"""
struct Polar{TF,TS<:Vector{TF}}
    lift::TF
    drag::TF
    pdrag::TF
    idrag::TF
    moment::TF
    surfacevelocity::TS
    surfacepressure::TS
end

# """
# ViscousSolution{}

# **Fields:**
#  - `panelgammas::Array{Float,2}` : \$\\gamma_0\$ and \$\\gamma_{90}\$ values at each airfoil node.
#  - `panelsources::Array{Float,2}` : source values at each airfoil node.
#  - `wakesources::Array{Float,2}` : source values at each wake node.
#  - `bodystrength::Array{Float}` : \$\\Psi_0\$ (constant stream function) 0 and 90 values.
# """
# struct ViscousSolution{TF,TD}
#     panelgammas::TF
#     panelsources::TF
#     wakesources::TF
#     bodystrength::TF
# end

#"""
#    Properties{TF}

#Thermodynamic properties for the viscous solution

#**Fields:**
# - `machinf::Float` : freestream mach number
# - `KTbeta::Float` : Karman-Tsien beta
# - `KTlambda::Float` : Karman-Tsien lambda
# - `H0::Float` : stagnation enthalpy
# - `sonic_cp::Float` : sonic cp
# - `rho0::Float` : stagnation density
# - `mu0::Float` : stagnation dynamic viscosity
#"""
#struct Properties{TF}
#    machinf::TF
#    KTbeta::TF
#    KTlambda::TF
#    H0::TF
#    sonic_cp::TF
#    rho0::TF
#    mu0::TF
#end

#"""
#    Parameters{TF}

#Solver Parameters.

#**Fields:**
# - `gamma_air::Float = 1.4` : ratio of specific heats for air
# - `eta_crit::Float = 9.0` : critical amplification factor
# - `eta_D::Float = 0.9` : wall/wake dissipation length ratio
# - `GA::Float = 6.7` : G - Beta locus A constant
# - `GB::Float = 0.75` : G - Beta locus B constant
# - `GC::Float = 18.0` : G - Beta locus C constant
# - `Klag::Float = 5.6` : shear lag constant
# - `Ctau::Float = 1.8` : shear stress initialization constant
# - `Etau::Float = 3.3` : shear stree initialization exponent
# - `rSu::Float = 0.35` : Sutherland temperature ratio
# - `fw::Float = 2.5` : wake gap continuation factor
# - `dw::Float = 1.0` : wake length, in airfoil chords
# - `epsilonw::Float = 1e-5` : first wake point offset, in airfoil chords
# - `iknowwhatimdoing::Bool` : boolean to silence warnings if you really know what you're doing.
# - `rhoinf::Float` : non-dimensional freestream density
# - `vinf::Float` : non-dimensional freestream velocity magnitude
# - `muinf::Float` : freestream dynamic viscosity
#"""
#struct Parameters{TF,TB}
#    gamma_air::TF
#    etacrit::TF
#    etaD::TF
#    GA::TF
#    GB::TF
#    GC::TF
#    Klag::TF
#    Ctau::TF
#    Etau::TF
#    rSu::TF
#    fw::TF
#    dw::TF
#    epsilonw::TF
#    iknowwhatimdoing::TB
#    rhoinf::TF
#    vinf::TF
#    muinf::TF
#end

#"""
#    defaultparameters(;kwargs)

#Initialized parameters to defaults, but allows selective user override through keyword arguments.

#**Keyword Arguments:**
# - `gamma_air::Float` : ratio of specific heats for air
# - `eta_crit::Float` : critical amplification factor
# - `eta_D::Float` : wall/wake dissipation length ratio
# - `GA::Float` : G - Beta locus A constant
# - `GB::Float` : G - Beta locus B constant
# - `GC::Float` : G - Beta locus C constant
# - `Klag::Float` : shear lag constant
# - `Ctau::Float` : shear stress initialization constant
# - `Etau::Float` : shear stree initialization exponent
# - `rSu::Float` : Sutherland temperature ratio
# - `fw::Float` : wake gap continuation factor
# - `dw::Float` : wake length, in airfoil chords
# - `epsilonw::Float` : first wake point offset, in airfoil chords
# - `iknowwhatimdoing::Bool` : boolean to silence warnings if you really know what you're doing.
# - `rhoinf::Float` : non-dimensional freestream density
# - `vinf::Float` : non-dimensional freestream velocity magnitude
# - `muinf::Float` : freestream dynamic viscosity
#"""
#function defaultparameters(;
#    gamma_air=1.4,
#    eta_crit=9.0,
#    eta_D=0.9,
#    GA=6.7,
#    GB=0.75,
#    GC=18.0,
#    Klag=5.6,
#    Ctau=1.8,
#    Etau=3.3,
#    rSu=0.35,
#    fw=2.5,
#    dw=1.0,
#    epsilonw=1e-5,
#    iknowwhatimdoing=false,
#    rhoinf=1.0,
#    vinf=1.0,
#    muinf=0.0,
#)

#    #check if vinf and rhoinf are not default
#    if !iknowwhatimdoing && (vinf != 1.0 || rhoinf != 1.0)
#        @warn(
#            "vinf and/or rhoinf have been changed from their defaults. Do you really want to do that??\n\nTo silence this warning in the future, set the field `iknowwhatimdoing=true`."
#        )
#    end

#    return Parameters(
#        gamma_air, eta_crit, eta_D, GA, GB, GC, Klag, Ctau, Etau, rSu, fw, dw, epsilonw
#    )
#end

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
    AxiSymPanel{TCP,TL,TNH,TB,TR}

Panel object for axisymmetric meshes.

**Fields:**
- `controlpoint::Array{Float}` : [x;r] coordinates of panel midpoint.
- `length::Float` : length of panel
- `normal::Array{Float}` : unit normal vector of panel (TODO: remove if unused)
- `beta::Float` : angle panel makes with positive x-axis (radians)
- `radiusofcurvature::Float` : the radius of curvature of the geometry at the panel control point. TODO: make sure this is actually correct with current implementation.
"""
struct AxiSymPanel{TCP,TL,TNH,TB,TR}
    controlpoint::TCP
    length::TL
    normal::TNH
    beta::TB
    radiusofcurvature::TR
end
