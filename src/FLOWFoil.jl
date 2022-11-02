module FLOWFoil

# - DEPENDENCIES
using LinearAlgebra
using FLOWMath
using SpecialFunctions

# - EXPORTS

#TYPES
# geometry types
export PlanarMesh, PlanarMeshSystem, AxiSymMesh
# problem types
export Problem, InviscidSolution
# output types
export Polar

#FUNCTIONS
# geometry functions
export generate_mesh, generate_axisym_mesh, position_coordinates, position_coordinates!
# inviscid solver functions
export solve, inviscid_post, calculate_stream_grid
# common airfoil parameterizations
export karman_trefftz, joukowsky, naca4

# - INCLUDED FILES

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
include("../common_parameterizations/bspline.jl") #REQUIRES UNREGISTERED PACKAGE

end
