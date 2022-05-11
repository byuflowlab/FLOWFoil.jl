#=
Compose.g.Type Definitions

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
    wake_nodes::Array{Array{TF,2}}
    wake_midpoints::Array{Array{TF,2}}
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

"""
    Freestream{TF}

Freestream Definition.

**Fields:**
 - 'reynolds::Vector{Float}' : Reynolds Numbers.
 - 'machs::Vector{Float}' : Mach Numbers.
 - 'anglesofattack::Vector{Float}' : Angles of attack.

"""
struct Freestream{TF}
    reynolds::Vector{TF}
    density::Vector{TF}
    dynamicviscosity::Vector{TF}
    mach::Vector{TF}
    anglesofattack::Vector{TF}
end

"""
    Problem{TM,TS,TB}

Problem definition and method selection.

**Fields:**
 - 'meshsystem::MeshSystem : Mesh System to solve.
 - 'panelshape::String' : Shape of panel to use, e.g., "le.g.r".
 - 'singularitytype::String' : Type of singularity to use, e.g., "source".
 - 'singularityorder::String' : Order of singularity distributions, e.g., "constant".
 - 'verbose::Bool' : Flag to print out verbose statements.
 - 'debug::Bool' : Flag to print out debugging statements.

"""
struct Problem{TM,TS,TB}
    meshsystem::TM
    panelshape::TS
    singularitytype::TS
    singularityorder::TS
    verbose::TB
    debug::TB
end

"""
    Solution{TM,TF}

Output object containing solution and useful items.

**Fields:**
 - 'meshsystem::MeshSystem' : Mesh System used in solution (potentially modified from input meshes).
 - 'strengthsvec::Array{Float,2}' : singularity strengths.
 - 'geocoeffmat::Array{Float,2}' : Geometric Coefficient Matrix used in solution.
 - 'bccoeffvec::Vector{Float}' : Boundary Coefficient Vector used in solution.
 - 'lift::Vector{Float}' : Lift Coefficients.
 - 'drag::Vector{Float}' : Total Drag Coefficients.
 - 'pdrag::Vector{Float}' : Pressure Drag Coefficients.
 - 'idrag::Vector{Float}' : Induced Drag Coefficients.
 - 'moment::Vector{Float}' : Moment Coefficients.

"""
struct Solution{TM,TF}
    meshsystem::TM
    strengthsvec::Array{TF}
    geocoeffmat::Array{TF,2}
    bccoeffvec::Vector{TF}
    lift::Vector{TF}
    drag::Vector{TF}
    pdrag::Vector{TF}
    idrag::Vector{TF}
    moment::Vector{TF}
end
