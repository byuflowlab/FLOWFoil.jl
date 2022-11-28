module FLOWFoil

# - DEPENDENCIES
using LinearAlgebra
using FLOWMath
using SpecialFunctions

# - EXPORTS

#TYPES
# geometry types
export PlanarMesh, PlanarMeshSystem, AxiSymMesh, AxiSymPanel
# problem types
export Problem, InviscidSolution
# output types
export PlanarPolar, AxiSymPolar

#FUNCTIONS
# geometry functions
export generate_mesh, generate_axisym_mesh, position_coordinates, position_coordinates!
# inviscid solver functions
export solve
# post processing functions
export get_planar_polar, get_axisymmetric_polar
# common airfoil parameterizations
export karman_trefftz, joukowsky, naca4#, gbs

# - INCLUDED FILES

# Geometry Generation and Modification
include("panel.jl")

# Singularity Distributions
include("singularity.jl")

# Inviscid Solver
include("inviscid_system.jl")

# Boundary Layer Integration
# include("viscous_system.jl")

# Solver
include("solve.jl")

# Post Processing
include("post_process.jl")

include("planar_post_process.jl")
include("axisymmetric_post_process.jl")

# Common Airfoil Parameterizations
include("airfoils/parameterizations/utils.jl")
include("airfoils/parameterizations/conformal_mapping.jl")
include("airfoils/parameterizations/naca.jl")
# include("airfoils/parameterizations/cst.jl")
# include("airfoils/parameterizations/parsec.jl")
# include("airfoils/parameterizations/bspline.jl") #REQUIRES UNREGISTERED PACKAGE

end
