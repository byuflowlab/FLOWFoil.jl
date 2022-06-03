#=
Compose Type Definitions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

#=
Organization TODO:

- User input = Problem type.
- problem includes airfoil coordinates and freestream definition (independant non-dim parameters)
- For each body, an airfoil mesh is created
- for each body, a wake mesh is created
- an inviscid system is used to pass inviscid coefficients and boundary conditions, etc.
- a viscous system is used to pass all the viscous coefficients, boundary layer numbers, etc.
- solution objects include everything needed to understand the flow field solved for and can be passed into post-process function to get all the xfoil-like outputs (also need convenience function to go from problem straight to cl, cd, cp, cm, etc.)
- debug option includes debug struct in solution.
- debug struct contains the system objects and whatever else could be useful.

=#

"""
    Problem{TF,TB}

Problem definition (geometry, operating point(s), and method selection) and output behavior.

**Fields:**
 - 'coordinates::Array{Float}' : x,y airfoil coordinates.
 - 'angleofattack::Float' : angle of attack to analyze.
 - 'reynolds::Float' : Reynolds number to analyze.
 - 'mach::Float' : Mach number to analyze.
 - 'viscous::Bool' : Flag to solve viscous or inviscid only
 - 'verbose::Bool' : Flag to print out verbose statements
 - 'debug::Bool' : Flag to save the system structs, etc.

"""
struct Problem{TM,TF,TB}
    meshes::Array{TM}
    angleofattack::TF
    reynolds::TF
    mach::TF
    viscous::TB
    verbose::TB
    debug::TB
end

function Problem(
    meshes, angleofattack, reynolds, mach=0.0; viscous=true, verbose=false, debug=false
)
    return Problem(meshes, angleofattack, reynolds, mach, viscous, verbose, debug)
end

"""
    BodyMesh{TF,TB}

Mesh for single body.

**Fields:**
 - 'airfoil_nodes::Array{Array{Float,2}}' : [x y] node (panel edge) locations for airfoil
 - 'chord::Float' : airfoil chord length
 - 'blunt_te::Bool' : boolean for whether or not the trailing edge is blunt or not.
**Assuptions:**
 - x and y coordinates start at the bottom trailing edge and proceed clockwise.

"""
struct BodyMesh{TF,TB}
    airfoil_nodes::Array{Matrix{TF}}
    chord::TF
    blunt_te::TB
    trailing_edge_gap::TF
    tdp::TF
    txp::TF
end

"""
    BodyMeshSystem{TF}

# System of meshes to solve.

# **Fields:**
#  - 'meshes::Array{Mesh}' : Array of mesh objects.
#  - 'scales::Vector{Float}' : Airfoil scaling factors.
#  - 'angles::Vector{Float}' : Airfoil angles of attack.
#  - 'locations::Array{Array{TF}}' : Array of leading edge locations.

"""
struct BodyMeshSystem{TM,TF}
    meshes::Vector{TM}
    scales::Vector{TF}
    angles::Vector{TF}
    locations::Vector{Vector{TF}}
end

"""
    WakeMesh{TF}

**Fields:**
 - 'wake_nodes::Array{Float,2}' : x,y wake panel node locations.
 - 'wake_midpoints::Array{Float,2}' : x,y wake panel center point locations.
"""
struct WakeMesh{TF}
    wake_nodes::Array{TF}
    wake_midpoints::Array{TF}
end

"""
    InviscidSystem{TF}

**Fields:**
 - 'vcoeffmat::Array{Float,2}' : Vortex Coefficient Matrix used in solution.
 - 'bccoeffvec::Array{Float,2}' : Boundary Coefficient Vector used in solution.
"""
struct InviscidSystem{TF,TI}
    vcoeffmat::Array{TF}
    bccoeffvec::Array{TF}
    Ns::Array{TI}
end

"""
**Fields:**
"""
struct ViscousSystem{TPa,TPr,TF}
    parameters::TPa
    properties::TPr
    vcoeffmat::Array{TF}
    scoeffmat::Array{TF}
    bccoeffvec::Array{TF}
end

"""
    InviscidSolution{TM,TF,TD}

**Fields:**
 - 'mesh::BodyMesh' : BodyMesh object describing airfoil nodes etc.
 - 'panelgammas::Array{Float,2}' : \$\\gamma_0\$ and \$\\gamma_{90}\$ values at each airfoil node.
 - 'psi0::Array{Float}' : \$\\Psi_0\$ (constant stream function) 0 and 90 values.
 - 'debug::Debug' : Debug object (or nothing) depending on debug flag in Problem object.
"""
struct InviscidSolution{TM,TF,TI,TD}
    meshes::Array{TM}
    panelgammas::Array{TF}
    psi0::Array{TF}
    Ns::Array{TI}
    debug::TD
end

"""
    ViscousSolution{}

**Fields:**
 - 'panelgammas::Array{Float,2}' : \$\\gamma_0\$ and \$\\gamma_{90}\$ values at each airfoil node.
 - 'panelsources::Array{Float,2}' : source values at each airfoil node.
 - 'wakesources::Array{Float,2}' : source values at each wake node.
 - 'psi0::Array{Float}' : \$\\Psi_0\$ (constant stream function) 0 and 90 values.

"""
struct ViscousSolution{TF,TD}
    panelgammas::Array{TF}
    panelsources::Array{TF}
    wakesources::Array{TF}
    psi0::Array{TF}
    debug::TD
end

# TODO: Need to figure out how to define/use debug object.
"""
    Debug{TIS,TVS}

**Fields:**
 - 'isystem::InviscidSystem' : Inviscid System Object.
 - 'vsystem::ViscousSystem' : Viscous System Object.
"""
struct Debug{TIS,TVS}
    isystem::TIS
    vsystem::TVS
end

"""
    Polar{TF}

**Fields:**
 - 'lift::Float' : Lift Coefficient.
 - 'drag::Float' : Total Drag Coefficient.
 - 'pdrag::Float' : Pressure Drag Coefficient.
 - 'idrag::Float' : Induced Drag Coefficient.
 - 'moment::Float' : Moment Coefficient.
 - 'surfacevelocity::Vector{Float}' : surface velocity distribution
 - 'surfacepressure::Vector{Float}' : surface pressure distribution
"""
struct Polar{TF}
    lift::TF
    drag::TF
    pdrag::TF
    idrag::TF
    moment::TF
    surfacevelocity::Vector{TF}
    surfacepressure::Vector{TF}
end

"""
    Properties{TF}

Thermodynamic properties for the viscous solution

**Fields:**
 - 'machinf::Float' : freestream mach number
 - 'KTbeta::Float' : Karman-Tsien beta
 - 'KTlambda::Float' : Karman-Tsien lambda
 - 'H0::Float' : stagnation enthalpy
 - 'sonic_cp::Float' : sonic cp
 - 'rho0::Float' : stagnation density
 - 'mu0::Float' : stagnation dynamic viscosity
"""
struct Properties{TF}
    machinf::TF
    KTbeta::TF
    KTlambda::TF
    H0::TF
    sonic_cp::TF
    rho0::TF
    mu0::TF
end

"""
    Parameters{TF}

Solver Parameters.

**Fields:**
 - 'gamma_air::Float = 1.4' : ratio of specific heats for air
 - 'eta_crit::Float = 9.0' : critical amplification factor
 - 'eta_D::Float = 0.9' : wall/wake dissipation length ratio
 - 'GA::Float = 6.7' : G - Beta locus A constant
 - 'GB::Float = 0.75' : G - Beta locus B constant
 - 'GC::Float = 18.0' : G - Beta locus C constant
 - 'Klag::Float = 5.6' : shear lag constant
 - 'Ctau::Float = 1.8' : shear stress initialization constant
 - 'Etau::Float = 3.3' : shear stree initialization exponent
 - 'rSu::Float = 0.35' : Sutherland temperature ratio
 - 'fw::Float = 2.5' : wake gap continuation factor
 - 'dw::Float = 1.0' : wake length, in airfoil chords
 - 'epsilonw::Float = 1e-5' : first wake point offset, in airfoil chords
 - 'iknowwhatimdoing::Bool' : boolean to silence warnings if you really know what you're doing.
 - 'rhoinf::Float' : non-dimensional freestream density
 - 'vinf::Float' : non-dimensional freestream velocity magnitude
 - 'muinf::Float' : freestream dynamic viscosity
"""
struct Parameters{TF,TB}
    gamma_air::TF
    etacrit::TF
    etaD::TF
    GA::TF
    GB::TF
    GC::TF
    Klag::TF
    Ctau::TF
    Etau::TF
    rSu::TF
    fw::TF
    dw::TF
    epsilonw::TF
    iknowwhatimdoing::TB
    rhoinf::TF
    vinf::TF
    muinf::TF
end

"""
    defaultparameters(;kwargs)

Initialized parameters to defaults, but allows selective user override through keyword arguments.

**Keyword Arguments:**
 - 'gamma_air::Float' : ratio of specific heats for air
 - 'eta_crit::Float' : critical amplification factor
 - 'eta_D::Float' : wall/wake dissipation length ratio
 - 'GA::Float' : G - Beta locus A constant
 - 'GB::Float' : G - Beta locus B constant
 - 'GC::Float' : G - Beta locus C constant
 - 'Klag::Float' : shear lag constant
 - 'Ctau::Float' : shear stress initialization constant
 - 'Etau::Float' : shear stree initialization exponent
 - 'rSu::Float' : Sutherland temperature ratio
 - 'fw::Float' : wake gap continuation factor
 - 'dw::Float' : wake length, in airfoil chords
 - 'epsilonw::Float' : first wake point offset, in airfoil chords
 - 'iknowwhatimdoing::Bool' : boolean to silence warnings if you really know what you're doing.
 - 'rhoinf::Float' : non-dimensional freestream density
 - 'vinf::Float' : non-dimensional freestream velocity magnitude
 - 'muinf::Float' : freestream dynamic viscosity
"""
function defaultparameters(;
    gamma_air=1.4,
    eta_crit=9.0,
    eta_D=0.9,
    GA=6.7,
    GB=0.75,
    GC=18.0,
    Klag=5.6,
    Ctau=1.8,
    Etau=3.3,
    rSu=0.35,
    fw=2.5,
    dw=1.0,
    epsilonw=1e-5,
    iknowwhatimdoing=false,
    rhoinf=1.0,
    vinf=1.0,
    muinf=0.0,
)

    #check if vinf and rhoinf are not default
    if !iknowwhatimdoing && (vinf != 1.0 || rhoinf != 1.0)
        @warn(
            "vinf and/or rhoinf have been changed from their defaults. Do you really want to do that??\n\nTo silence this warning in the future, set the field 'iknowwhatimdoing=true'."
        )
    end

    return Parameters(
        gamma_air, eta_crit, eta_D, GA, GB, GC, Klag, Ctau, Etau, rSu, fw, dw, epsilonw
    )
end
