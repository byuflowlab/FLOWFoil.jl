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
    Problem{TF}

Problem definition (geometry, operating point(s), and method selection) and output behavior.

**Fields:**
 - `mesh::Array{Mesh}` : Array of mesh objects
 - `aoa::Float` : angle of attack to analyze.
 - `reynolds::Float` : Reynolds number to analyze.
 - `mach::Float` : Mach number to analyze.
 - `viscous::Bool` : Flag to solve viscous or inviscid only
- `solver::ProblemType` : Analysis method to use for problem.
"""
struct Problem{TF}
    nbody::Int
    angle_of_attack::Vector{TF}
    reynolds::Vector{TF}
    mach::Vector{TF}
    viscous::Bool
    type::ProblemType
end

"""
    define_problem(problemtype::ProblemType, coordinates, angle_of_attack, reynolds, mach)

Defines a Problem object to be used throughout setup, solution, and post-processing.

**Arguments:**
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)
- `angle_of_attack::Vector{Float}` : Vector of angles of attack (may be a single float as well)
- `reynolds::Vector{Float}` : Vector of reynolds numbers (may be a single float as well)
- `mach::Vector{Float}` : Vector of mach numbers (may be a single float as well)
- `problemtype::ProblemType` : Type of problem to solve (planar, axisymmetric, or periodic).

**Returns:**
- `problem::Problem` : Problem object
"""
function define_problem(
    problemtype::ProblemType, coordinates, angle_of_attack, reynolds, mach
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
    if length(angle_of_attack) == 1
        aoa = [angle_of_attack[1]]
    else
        aoa = angle_of_attack
    end

    # - Make Reynolds Number a Vector (if not one already) - #
    if reynolds == nothing
        reynolds_numbers = [reynolds[1]]
    elseif length(reynolds) == 1
        reynolds_numbers = [reynolds[1]]
    else
        reynolds_nubmers = reynolds
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

# - Simple Implementation (Like XFoil) - #
function define_problem(coordinates, angle_of_attack, reynolds=nothing, mach=nothing;)

    # If no solver type, use planar
    return define_problem(
        Planar(), coordinates, angle_of_attack; reynolds=reynolds, mach=mach
    )
end
