
######################################################################
#                                                                    #
#                         Singularity Types                          #
#                                                                    #
######################################################################

abstract type Singularity end
abstract type Order end

### --- Singularity Types --- ###

struct Source{TO} <: Singularity where {TO<:Order}
    order::TO
end

struct Vortex{TO} <: Singularity where {TO<:Order}
    order::TO
end

struct Doublet{TO} <: Singularity where {TO<:Order}
    order::TO
end

### --- Singularity Order Types --- ###

struct Constant <: Order end

struct Linear <: Order end

struct Quadratic <: Order end

struct Spline <: Order end

######################################################################
#                                                                    #
#                       BOUNDARY CONDITIONS                          #
#                                                                    #
######################################################################

abstract type BoundaryCondition end

# Neumann BC defines velocity, viz., the derivative of the potential
# (this is the no flow through condition)
struct Neumann <: BoundaryCondition end

# Dirichlet BC defines potential (zero potential inside the body)
struct Dirichlet <: BoundaryCondition end

# Robin BC is a weighted combination of Dirichlet and Neumann over the entire boundary
# Will take some extra stuff to get working
struct Robin <: BoundaryCondition end

# Mixed BC is where Neumann and Dirichlet BCs are used at disparate locations on the boundary
# Also will take some extra stuff to get working
struct Mixed <: BoundaryCondition end

######################################################################
#                                                                    #
#                          PROBLEM TYPES                             #
#                                                                    #
######################################################################

# Define abstract type of which each solver struct will be a subtype
abstract type ProblemType end

# - Planar (2D) - #
# This is the standard airfoil analysis type (similar to XFoil)
"""
    Planar()

One of the available sub-types of ProblemType.  Indicates use of Planar (2D) analysis method.
"""
struct PlanarProblem{TS,TBC} <: ProblemType where {TS<:Singularity,TBC<:BoundaryCondition}
    singularity::TS
    boundary::TBC
end
# NOTE: Currently, the only Planer implementation is the mfoil solver.
# TODO: eventually, allow for user defined panel and singularity types

# - Axisymmetric - #
# The Axisymmetric solver type is used for bodies of revolution and annular airfoils (ducts)
"""
    Axisymmetric(body_of_revolution)

One of the available sub-types of ProblemType.  Indicates use of Axisymmetric analysis method.

**Fields:**
- `body_of_revolution::Vector{Bool}` : Array of bools indicating whether the associated body is a body of revolution or not. (If not, it will be considered an annular airfoil, i.e., a duct.)
"""
struct AxisymmetricProblem <: ProblemType
    body_of_revolution::Vector{Bool}
end

# - Periodic (Cascade) - #
# The periodic type is specifically used for cascade analysis.
"""
    Periodic()

One of the available sub-types of ProblemType.  Indicates use of Periodic (cascade) analysis method.
"""
struct PeriodicProblem <: ProblemType end
