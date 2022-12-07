#=

Abstract types and structs used in dispatching.

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                         Singularity Types                          #
#                                                                    #
######################################################################

abstract type Singularity end
abstract type Order end

### --- Singularity Types --- ###

"""
    Source <: Singularity

One of the types of singularities available.

**Fields:**
- `order::Order` : order of singularities
"""
struct Source{TO} <: Singularity where {TO<:Order}
    order::TO
end

"""
    Vortex <: Singularity

One of the types of singularities available.

**Fields:**
- `order::Order` : order of singularities
"""
struct Vortex{TO} <: Singularity where {TO<:Order}
    order::TO
end

"""
    Doublet <: Singularity

One of the types of singularities available.

**Fields:**
- `order::Order` : order of singularities
"""
struct Doublet{TO} <: Singularity where {TO<:Order}
    order::TO
end

### --- Singularity Order Types --- ###

"""
    Constant <: Order

One of the singularity orders available.
"""
struct Constant <: Order end

"""
    Linear <: Order

One of the singularity orders available.
"""
struct Linear <: Order end

"""
    Quadratic <: Order

One of the singularity orders available.
"""
struct Quadratic <: Order end

"""
    Spline <: Order

One of the singularity orders available.
Used for **I**so**g**eometric **A**nalyis (IgA) implementations.
"""
struct Spline <: Order end

######################################################################
#                                                                    #
#                       BOUNDARY CONDITIONS                          #
#                                                                    #
######################################################################

abstract type BoundaryCondition end

"""
    Neumann <: BoundaryCondition

One of the available boundary conditions.

The Neumann boundary condition defines velocity, viz., the derivative of the potential, on the boundary. (This is the no flow through condition.)
"""
struct Neumann <: BoundaryCondition end

"""
    Dirichlet <: BoundaryCondition

One of the available boundary conditions.

The Dirichlet boundary condition defines the potential on the boundary. (This is the zero potential inside the body condition.)
"""
struct Dirichlet <: BoundaryCondition end

"""
    Robin

One of the available boundary conditions.

The Robin boundary condition is a weighted combination of Dirichlet and Neumann conditions over the entire boundary.
"""
struct Robin <: BoundaryCondition end
# NOTE: will need to add some fields here if used in the future to contain weighting information

"""
    Mixed

One of the available boundary conditions.

The Mixed boundary condition is where a combination of Neumann and Dirichlet conditions are used at disparate locations along the boundary.
"""
struct Mixed <: BoundaryCondition end
# NOTE: will need to add some fields here if used in the future to contain location information

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

**Fields:**
- `singularity::Singularity` : The type of singularities to be used in the problem.
- `boundary::BoundaryCondition` : The type of boundary condition to be used in the problem.
"""
struct PlanarProblem{TS,TBC} <: ProblemType where {TS<:Singularity,TBC<:BoundaryCondition}
    singularity::TS
    boundary::TBC
end
# NOTE: Currently, the only Planar implementation is the mfoil solver.
# NOTE: Eventually, will likely need to update types here to be vectors to allow for combinations of singularities (or different singularities on different bodies).  (boundary conditions already have mixed types)

# - Axisymmetric - #
# The Axisymmetric solver type is used for bodies of revolution and annular airfoils (ducts)
"""
    Axisymmetric(body_of_revolution)

One of the available sub-types of ProblemType.  Indicates use of Axisymmetric analysis method.

**Fields:**
- `body_of_revolution::Vector{Bool}` : Array of bools indicating whether the associated body is a body of revolution or not. (If not, it will be considered an annular airfoil, i.e., a duct.)
"""
struct AxisymmetricProblem{TS,TBC} <:
       ProblemType where {TS<:Singularity,TBC<:BoundaryCondition}
    body_of_revolution::Vector{Bool}
    singularity::TS
    boundary::TBC
end

# - Periodic (Cascade) - #
# The periodic type is specifically used for cascade analysis.
"""
    Periodic()

One of the available sub-types of ProblemType.  Indicates use of Periodic (cascade) analysis method.
"""
struct PeriodicProblem <: ProblemType end
