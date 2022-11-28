module FLOWFoil

##### ----- DEPENDENCIES ----- #####
using LinearAlgebra
using FLOWMath
using SpecialFunctions

##### ----- EXPORTS ----- #####

### --- TYPES --- ###

# - Geometry - #
export PlanarMesh, PlanarMeshSystem, AxiSymMesh, AxiSymPanel

# - Problem - #
export Problem, InviscidSolution

# - Polar - #
export PlanarPolar, AxiSymPolar

### --- FUNCTIONS --- ###

# - Geometry - #
export generate_mesh, generate_axisym_mesh, position_coordinates, position_coordinates!

# - Solution - #
export solve

# - Polar - #
export get_planar_polar, get_axisymmetric_polar

# - Airfoils - #
export karman_trefftz, joukowsky, naca4#, gbs

##### ----- INCLUDES ----- #####

# Problem Object Definition
include("problem.jl")

# Panel Geometry Functions
include("panel.jl")

# Singularity Functions
include("singularity.jl")

# Linear System Generation
include("system.jl")

# Linear System Solve
include("solve.jl")

# Solution Post Processing
include("post_process.jl")

### --- Airfoil Parameterizations and Manipulations --- ###

# - MANIPULATIONS - #
include("airfoils/parameterizations/utils.jl")

# - PARAMETERIZAIONS - #

# Conformal Mapping
include("airfoils/parameterizations/conformal_mapping.jl")

# NACA 4-Series
include("airfoils/parameterizations/naca.jl")

# Class-Shape Transformation
# include("airfoils/parameterizations/cst.jl")

# PARSEC
# include("airfoils/parameterizations/parsec.jl")

# B-Spline
# include("airfoils/parameterizations/bspline.jl") #REQUIRES UNREGISTERED PACKAGE

end
