#=

Problem Definition Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                          PROBLEM TYPES                             #
#                                                                    #
######################################################################

# Define abstract type of which each solver struct will be a subtype
abstract type ProblemType end

# - Planar (2D) - #
# This is the standard airfoil analysis type (similar to XFoil)
struct Planar <: ProblemType end
# NOTE: Currently, the only Planer implementation is the mfoil solver.
# TODO: eventually, allow for user defined panel and singularity types

# - Axisymmetric - #
# The Axisymmetric solver type is used for bodies of revolution and annular airfoils (ducts)
struct Axisymmetric{TB} <: ProblemType
    body_of_refolution::TB
end
# QUESTION FOR TAYLOR: Here TB could be a single bool or an array of bools.  Should I keep it as TB, or should I always make it an array of bools explicitly?

# - Periodic (Cascade) - #
# The periodic type is specifically used for cascade analysis.
struct Periodic <: ProblemType end

######################################################################
#                                                                    #
#                         PROBLEM DEFINITION                         #
#                                                                    #
######################################################################

# - Type Definition - #
"""
    Problem{TM,TF,TB}

Problem definition (geometry, operating point(s), and method selection) and output behavior.

**Fields:**
 - `mesh::Array{Mesh}` : Array of mesh objects
 - `aoa::Float` : angle of attack to analyze.
 - `reynolds::Float` : Reynolds number to analyze.
 - `mach::Float` : Mach number to analyze.
 - `viscous::Bool` : Flag to solve viscous or inviscid only
- `solver::ProblemType` : Analysis method to use for problem.
"""
struct Problem{TM,TF,TB}
    mesh::TM
    aoa::TF
    reynolds::TF
    mach::TF
    viscous::TB
    type::ProblemType
end

# - Basic Function - #
"""
"""
function define_problem(
    ::ProblemType,
    coordinates,
    angle_of_attack,
    reynolds,
    mach;
    scaling=[],
    translation=[],
    flip=[],
) end

# - Simple Implementation (Like XFoil) - #
function define_problem(
    coordinates,
    angle_of_attack,
    reynolds=nothing,
    mach=nothing;
    scaling=[],
    translation=[],
    flip=[],
)

    # If no solver type, use planar
    return define_problem(
        Planar(),
        coordinates,
        angle_of_attack;
        reynolds=reynolds,
        mach=mach;
        scaling=scaling,
        translation=translation,
        flip=flip,
    )
end

# - Planar Implementation - #
"""
"""
function define_problem(
    ::Planar;
    coordinates,
    angle_of_attack,
    reynolds,
    mach;
    scaling=[],
    translation=[],
    flip=[],
) end

# - Axisymmetric Implementation - #
"""
"""
function define_problem(
    ::Axisymmetric;
    coordinates,
    angle_of_attack,
    reynolds,
    mach;
    scaling=[],
    translation=[],
    flip=[],
) end

# - Periodic Implementation - #
"""
"""
function define_problem(
    ::Periodic;
    coordinates,
    angle_of_attack,
    reynolds,
    mach;
    scaling=[],
    translation=[],
    flip=[],
) end
