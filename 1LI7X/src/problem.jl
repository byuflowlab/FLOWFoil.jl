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

"""
    define_problem(problemtype::ProblemType, coordinates, flow_angle, reynolds, mach)

Defines a Problem object to be used for setup and post-processing.

**Arguments:**
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `flow_angle::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)
- `problemtype::ProblemType` : Type of problem to solve (planar, axisymmetric, or periodic).

**Returns:**
- `problem::Problem` : Problem object
"""
function define_problem(
    problemtype::ProblemType, coordinates, flow_angle, reynolds, mach
)

    # - Get Number of Bodies - #
    #= NOTE:
        The size function in the first if will use the overloaded version that returns the length of the tuple if coordinates is of type Tuple.  See overloaded function at the end of FLOWFoil.jl
    =#
    if length(size(coordinates)) == 1
        nbody = length(coordinates)
    else
        nbody = size(coordinates)[2]
    end

    # - Make Angles of Attack a Vector (if not one already) - #
    if length(flow_angle) == 1
        aoa = [flow_angle[1]]
    else
        aoa = flow_angle
    end

    # - Make Reynolds Number a Vector (if not one already) - #
    if reynolds == nothing
        reynolds_numbers = [reynolds[1]]
    elseif length(reynolds) == 1
        reynolds_numbers = [reynolds[1]]
    else
        reynolds_numbers = reynolds
    end

    # - Make Mach Number a Vector (if not one already) - #
    if mach == nothing
        mach_numbers = [mach[1]]
    elseif length(mach) == 1
        mach_numbers = [mach[1]]
    else
        mach_numbers = mach
    end

    # - Must be Viscous if Reynolds Numbers are Greater Than Zero - #
    if all(x -> x != nothing, reynolds_numbers) && all(x -> x > 0.0, reynolds_numbers)
        viscous = true
    elseif all(x -> x == nothing, reynolds_numbers) || all(x -> x <= 0.0, reynolds_numbers)
        viscous = false
    else
        @error "Cannot mix viscid and inviscid analyses. For viscous analsis sets, all reynolds numbers must be greater than 0.0. For inviscid analysis sets, it is best to set reynolds to -1.0."
    end

    # - RETURN - #
    return Problem(nbody, aoa, reynolds_numbers, mach_numbers, viscous, problemtype)
end
