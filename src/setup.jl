#=

Problem Definition Functions

Authors: Judd Mehr,

=#

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
- `solver::AnalysisType` : Analysis method to use for problem.
"""
struct Problem{TM,TF,TB}
    mesh::TM
    aoa::TF
    reynolds::TF
    mach::TF
    viscous::TB
    solver::AnalysisType
end

# - Basic Function - #
"""
"""
function define_problem(
    ::AnalysisType,
    coordinates,
    angle_of_attack,
    reynolds,
    mach;
    scaling=[],
    translation=[],
    flip=[],
) end

# - Simple Implementation (Like XFoil) - #
"""
"""
function define_problem(
    coordinates,
    angle_of_attack,
    reynolds=-1.0,
    mach=-1.0;
    scaling=[],
    translation=[],
    flip=[],
)

    # If no solver type, use planar
    return define_problem(
        ::Planar,
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
