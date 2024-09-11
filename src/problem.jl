#=

Problem Definition Types and Functions

"Problems" allow users to choose what angles of attack, Reynolds numbers, and Mach numbers they want to analyze their airfoil system for.
Users also define what type of problem they are solve (e.g. planar, axisymmetric, or periodic).

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                         PROBLEM DEFINITION                         #
#                                                                    #
######################################################################

# - Type Definition - #
"""
    Problem{TF}

Problem definition and output behavior.

**Fields:**
- `nbody::Int` : Number of bodies in system to analyze.
- `flow_angle::Vector{Float}` : angle(s) of attack to analyze.
- `reynolds::Vector{Float}` : Reynolds number(s) to analyze.
- `mach::Vector{Float}` : Mach number(s) to analyze.
- `viscous::Bool` : Flag to solve viscous or inviscid only
- `method::ProblemType` : Analysis method to use for problem.
"""
struct Problem{TF}
    nbody::Int
    flow_angle::Vector{TF}
    reynolds::Vector{TF}
    mach::Vector{TF}
    viscous::Bool
    method::ProblemType
end

