module FLOWFoil

# DEPENDENCIES
using LinearAlgebra

# EXPORTS

#TYPES
# geometry types
export Mesh, MeshSystem
# problem types
export Problem, Solution

#FUNCTIONS
# geometry functions
export generatemesh, position_meshes!
# inviscid solver functions
export solve, inviscid_post
# common airfoil parameterizations
export karman_trefftz, joukowsky, naca4

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
# include("viscous_system.jl")

# Solver
include("solve.jl")

# Post Processing
include("post_processes.jl")

# Common Airfoil Parameterizations
include("../common_parameterizations/convenience_functions.jl")
include("../common_parameterizations/conformal_mapping.jl")
include("../common_parameterizations/naca.jl")

end
