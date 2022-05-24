#=
Compose Type Definitions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    Mesh{TF}

Mesh for single body.

**Fields:**
 - 'airfoil_nodes::Array{Array{Float,2}}' : [x y] node (panel edge) locations for airfoil
 - 'wake_nodes::Array{Array{Float,2}}' : [x y] node (panel edge) locations for wake
 - 'wake_midpoints::Array{Array{Float,2}}' : [x y] wake panel midpoint locations
 - 'blunt_te::Bool' : boolean for whether or not the trailing edge is blunt or not.
**Assuptions:**
 - x and y coordinates start at the bottom trailing edge and proceed clockwise.
 - x and y coordinates are normalized such that the airfoil has chord length 1.0.

"""
struct Mesh{TF}
    airfoil_nodes::Array{Array{TF,2}}
    #    wake_nodes::Array{Array{TF,2}}
    #   wake_midpoints::Array{Array{TF,2}}
    blunt_te::Bool
end

"""
    MeshSystem{TF}

System of meshes to solve.

**Fields:**
 - 'meshes::Array{Mesh}' : Array of mesh objects.
 - 'scales::Vector{Float}' : Airfoil scaling factors.
 - 'angles::Vector{Float}' : Airfoil angles of attack.
 - 'locations::Array{Array{TF}}' : Array of leading edge locations.

"""
struct MeshSystem{TM,TF}
    meshes::Vector{TM}
    scales::Vector{TF}
    angles::Vector{TF}
    locations::Vector{Vector{TF}}
end

# TODO: probably make these type unions and set things up to loop through fields that are vectors rather than floats
"""
    Freestream{TF}

Freestream Definition.

**Fields:**
 - 'reynolds::Vector{Float}' : Reynolds Numbers
 - 'density::Vector{Float}' : air density
 - 'dynamicviscosity::Vector{Float}' : air dynamic viscosity
 - 'mach::Vector{Float}' : Mach Numbers
 - 'angleofattack::Vector{Float}' : Angles of attack in degrees

"""
struct Freestream{TF}
    reynolds::Vector{TF}
    density::Vector{TF}
    dynamicviscosity::Vector{TF}
    mach::Vector{TF}
    angleofattack::Vector{TF}
end

"""
    Problem{TM,TS,TB}

Problem definition and method selection.

**Fields:**
 - 'meshsystem::MeshSystem' : Mesh System to solve
 - 'freestream::FreeStream' : Freestream parameters
 - 'viscous::Bool' : Flag to solve viscous or inviscid only
 - 'verbose::Bool' : Flag to print out verbose statements
 - 'debug::Bool' : Flag to print out debugging statements

"""
struct Problem{TM,TS,TB}
    meshsystem::TM
    freestream::TS
    viscous::TB
    verbose::TB
    debug::TB
end

# TODO: Need an inviscid solution object as well as a viscous solution object.  They will have some similar fields and some different.
# # Also probably want to just have normal stuff in the outputs rather than every possible thing.  Probably create another type to house all of the minutia and if debug flag is true, also update and output that detailed object (come up with a reasonable name other than "Innards"
"""
    Solution{TM,TF}

Output object containing solution and useful items.

**Fields:**
 - 'meshsystem::MeshSystem' : Mesh System used in solution (potentially modified from input meshes).
 - 'strengthsvec::Array{Float,2}' : singularity strengths.
 - 'vcoeffmat::Array{Float,2}' : Vortex Coefficient Matrix used in solution.
 - 'scoeffmat::Array{Float,2}' : Source Coefficient Matrix used in solution.
 - 'bccoeffvec::Vector{Float}' : Boundary Coefficient Vector used in solution.
 - 'lift::Vector{Float}' : Lift Coefficients.
 - 'drag::Vector{Float}' : Total Drag Coefficients.
 - 'pdrag::Vector{Float}' : Pressure Drag Coefficients.
 - 'idrag::Vector{Float}' : Induced Drag Coefficients.
 - 'moment::Vector{Float}' : Moment Coefficients.
 - 'surfacevelocity::Vector{Float}' : surface velocity distribution
 - 'surfacepressure::Vector{Float}' : surface pressure distribution
"""
struct Solution{TM,TF}
    meshsystem::TM
    strengthsvec::Array{TF}
    vcoeffmat::Array{TF,2}
    scoeffmat::Array{TF,2}
    bccoeffvec::Vector{TF}
    lift::Vector{TF}
    drag::Vector{TF}
    pdrag::Vector{TF}
    idrag::Vector{TF}
    moment::Vector{TF}
    surfacevelocity::Vector{TF}
    surfacepressure::Vector{TF}
end

"""
    Parameters{TF}

 - `gamma_air::Float' : ratio of specific heats for air
 - `eta_crit::Float' : critical amplification factor
 - `eta_D::Float' : wall/wake dissipation length ratio
 - `GA::Float' : G - Beta locus A constant
 - `GB::Float' : G - Beta locus B constant
 - `GC::Float' : G - Beta locus C constant
 - `Klag::Float' : shear lag constant
 - `Ctau::Float' : shear stress initialization constant
 - `Etau::Float' : shear stree initialization exponent
 - `rSu::Float' : Sutherland temperature ratio
 - `fw::Float' : wake gap continuation factor
 - `dw::Float' : wake length, in airfoil chords
 - `epsilonw::Float' : first wake point offset, in airfoil chords
"""
struct Parameters{TF}
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
end

"""
    defaultparameters()

Initializes Parameters struct with defaults, see `Parameters` docstring.
"""
function defaultparameters()
    return Parameters(1.4, 9.0, 0.9, 6.7, 0.75, 18.0, 5.6, 1.8, 3.3, 0.35, 2.5, 1.0, 1e-5)
end
