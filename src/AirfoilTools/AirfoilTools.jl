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

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

# - Parameterizations - #
export NACA4, naca4, determine_naca4
export BasicBSpline, basic_bspline#, determine_basic_bspline
export KarmanTrefftz, karman_trefftz#, determine_karman_trefftz
export Joukowsky, joukowsky#, determine_joukowsky
export CST, cst, determine_cst
export PARSEC, parsec, determine_parsec
export PARSECStandard, parsec_standard, determine_parsec_standard

# - Manipulations - #
export whole_cosine_spacing, split_cosine_spacing, repanel_airfoil, refine_trailing_edge
export split_upper_lower
export flip!, zero_z_te!, rotate_coordinates!, normalize_coordinates!, position_coordinates!

#---------------------------------#
#            INCLUDES             #
#---------------------------------#

##### ----- PARAMETERIZATIONS ----- #####

# NACA 4-series
include("parameterizations/naca.jl")

# Class Shape Transformation
include("parameterizations/cst.jl")

# Generalized B-Spline (Rajnarayan)
include("parameterizations/basic_bspline.jl")

# Parametric Section
include("parameterizations/parsec.jl")

# Conformal Mapping
include("parameterizations/conformal_mapping.jl")

##### ----- GEOMTERY MANIPULATIONS ----- #####

# Splits, etc.
include("geometry_manipulations/deconstruction.jl")

# Rotations, Scaling, Transformations, etc.
include("geometry_manipulations/transformation.jl")

# Obtain additional information about coordinates
include("geometry_manipulations/inspection.jl")

# Interpolation, Smoothing, Reversing, Refinement, etc.
include("geometry_manipulations/redefinition.jl")

end
