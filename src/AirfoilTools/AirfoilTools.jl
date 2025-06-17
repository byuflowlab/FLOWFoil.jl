module AirfoilTools
#=

Collection of various airfoil geometry generation and manipulation functions.

Authors: Judd Mehr, Taylor McDonnell, Adam Cardoza, Justin Hawkins,

=#

#---------------------------------#
#           ABSTRACT TYPE         #
#---------------------------------#

abstract type AirfoilGeometry end

#---------------------------------#
#          DEPENDENCIES           #
#---------------------------------#
using LsqFit
using FLOWMath
using NURBS
using StaticArrays
using NLsolve
import ImplicitAD as iad

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

# - Parameterizations - #
export NACA4, naca4, determine_naca4
export naca65, naca65_scaled
export BasicBSpline, basic_bspline#, determine_basic_bspline
export KarmanTrefftz, karman_trefftz#, determine_karman_trefftz
export Joukowsky, joukowsky, joukowsky_flow#, determine_joukowsky
export CST, cst, determine_cst
export CircularArcCST, circular_arc_cst
export PARSEC, parsec, determine_parsec
export ModifiedPARSEC, modified_parsec, determine_modified_parsec

# - Manipulations - #
export whole_cosine_spacing,
    split_cosine_spacing, repanel_airfoil, refine_trailing_edge, linear_transform
export split_upper_lower
export flip!, zero_y_te!, rotate_coordinates!, normalize_coordinates!, position_coordinates!

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- PARAMETERIZATIONS ----- #####

# NACA 4-series
include("parameterizations/naca4.jl")

# NACA 65-series
include("parameterizations/naca65.jl")

# Class Shape Transformation (CST)
include("parameterizations/cst.jl")

# Circular Arc Camber line CST
include("parameterizations/circular_arc_cst.jl")

# Generalized B-Spline (Rajnarayan)
include("parameterizations/basic_bspline.jl")

# Parametric Section
include("parameterizations/parsec.jl")

# Joukowsky
include("parameterizations/joukowsky.jl")

# Karman_Trefftz
include("parameterizations/karman_trefftz.jl")

##### ----- GEOMTERY MANIPULATIONS ----- #####

# Splits, etc.
include("geometry_manipulations/deconstruction.jl")

# Rotations, Scaling, Transformations, etc.
include("geometry_manipulations/transformation.jl")

# Obtain additional information about coordinates
include("geometry_manipulations/inspection.jl")

# Interpolation, Smoothing, Reversing, Refinement, etc.
include("geometry_manipulations/redefinition.jl")

# Utilities
include("geometry_manipulations/utils.jl")

end
