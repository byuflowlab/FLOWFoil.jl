module FLOWFoil

# DEPENDENCIES
using LinearAlgebra

# EXPORTS

#TYPES
# geometry types
export Mesh, MeshSystem
# operating point types
export Freestream
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

# Inviscid Solver
include("inviscid_analysis.jl")

# Boundary Layer Integration
include("viscous_analysis.jl")

# Common Airfoil Parameterizations
include("../common_parameterizations/conformal_mapping.jl")

end
