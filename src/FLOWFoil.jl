module FLOWFoil

#---------------------------------#
#           DEPENDENCIES          #
#---------------------------------#
using LinearAlgebra
using FLOWMath
using SpecialFunctions

#---------------------------------#
#             INCLUDES            #
#---------------------------------#

##### ----- CORE FUNCTIONALITY ----- #####

# Dispatch Types
include("dispatch_types.jl")

# Problem Object Definition
include("problem.jl")

# Panel Definition
include("panel.jl")

# Mesh Generation
include("mesh.jl")

# Geometry Utilities
include("geometry.jl")

# Singularity Functions
include("singularity.jl")

# Linear System Generation
include("system.jl")

# Linear System Solve
include("solve.jl")

# # Solution Post Processing
# include("post_process.jl")

# Convenience Functions
include("convenience_functions.jl")

##### ----- AIRFOIL PARAMETERIZATIONS AND MANIPULATIONS ----- #####

# - MANIPULATIONS - #
include("airfoils/parameterizations/utils.jl")

# - PARAMETERIZAIONS - #

# Conformal Mapping
# include("airfoils/parameterizations/conformal_mapping.jl")

# NACA 4-Series
include("airfoils/parameterizations/naca.jl")

# Class-Shape Transformation
# include("airfoils/parameterizations/cst.jl")

# PARSEC
# include("airfoils/parameterizations/parsec.jl")

# B-Spline
# include("airfoils/parameterizations/bspline.jl") #REQUIRES UNREGISTERED PACKAGE

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

### --- TYPES --- ###
# Singularities
export Source, Doublet, Vortex
export Constant, Linear, Quadratic, Spline

# Boundary Conditions
export Neumann, Dirichlet, Robin, Mixed

# Problem
export Problem, PlanarProblem, AxisymmetricProblem, PeriodicProblem

# Panels
export PlanarFlatPanel, AxisymmetricFlatPanel

# Mesh
export generate_mesh

# System
export generate_inviscid_system

### --- FUNCTIONS --- ###

# Convenience Functions
export solve

# Problem
export define_problem

# Panels
export generate_panels

# Airfoil Parameterizations
export naca4

#---------------------------------#
#       OVERLOADED FUNCTIONS      #
#---------------------------------#

import Base.size
function size(t::Tuple)
    return length(t)
end

end
