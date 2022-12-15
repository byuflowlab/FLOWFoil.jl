module FLOWFoil

#---------------------------------#
#           DEPENDENCIES          #
#---------------------------------#
using LinearAlgebra
using FLOWMath
using ImplicitAD
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

# Solution Post Processing
include("post_process.jl")

# Convenience Functions
include("convenience_functions.jl")

# Utility Functions
include("utils.jl")

##### ----- AIRFOIL PARAMETERIZATIONS AND MANIPULATIONS ----- #####

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

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

##### ----- TYPES ----- #####

# Singularities
export Source, Doublet, Vortex
export Constant, Linear, Quadratic, Spline

# Boundary Conditions
export Neumann, Dirichlet, Robin, Mixed

# Problem
export Problem, PlanarProblem, AxisymmetricProblem, PeriodicProblem

# Panels
export PlanarFlatPanel, AxisymmetricFlatPanel

# System
export InviscidSystem

##### ----- FUNCTIONS ----- #####

# Convenience Functions
export solve

# Problem
export define_problem

# Panels
export generate_panels

# Mesh
export generate_mesh

# System
export generate_inviscid_system

# Post Process
export post_process

# Airfoil Parameterizations
export naca4

#---------------------------------#
#       OVERLOADED FUNCTIONS      #
#---------------------------------#

#= NOTE:
    Used in problem definition function to help count number of bodies, the coordinates of which are a tuple of vectors if multiple bodies are being analyzed together.
=#
import Base.size
function size(t::Tuple)
    return length(t)
end

end
