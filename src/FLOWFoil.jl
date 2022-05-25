module FLOWFoil

# DEPENDENCIES
using LinearAlgebra

# EXPORTS

#TYPES
# geometry types
export Mesh, MeshSystem
# operating point types
export Freestream, Parameters
# problem types
export Problem, Solution

#FUNCTIONS
# geometry functions
export generatemesh

# common airfoil parameterizations
export karmantrefftz, joukowsky

# INCLUDED FILES

# Types
include("types.jl")

# Geometry Generation and Modification
include("geometry.jl")

# Singularity Distributions
include("singularities.jl")

# Inviscid Solver
include("inviscid_system.jl")

# Boundary Layer Integration
include("viscous_system.jl")

# Solver
include("solve.jl")

# Common Airfoil Parameterizations
include("../common_parameterizations/conformal_mapping.jl")
include("../common_parameterizations/naca.jl")

end
