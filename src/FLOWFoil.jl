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

##### ----- AirfoilTools ----- #####

include("AirfoilTools/AirfoilTools.jl")
const at = AirfoilTools
export AirfoilTools

##### ----- Methods ----- #####

abstract type Method end

# Mfoil (Xfoil)
include("mfoil/method.jl")
include("mfoil/geometry_utils.jl")
include("mfoil/panel_geometry.jl")
include("mfoil/system_geometry.jl")
include("mfoil/singularities.jl")
include("mfoil/system_matrices.jl")
include("mfoil/solve.jl")
include("mfoil/post_process.jl")

# Lewis (Axisymmetric)
include("lewis/method.jl")
include("lewis/geometry_utils.jl")
include("lewis/panel_geometry.jl")
include("lewis/system_geometry.jl")
include("lewis/singularities.jl")
include("lewis/system_matrices.jl")
include("lewis/solve.jl")
include("lewis/post_process.jl")

# Martensen (Planar AND/OR Periodic)
include("martensen/method.jl")
include("martensen/panel_geometry.jl")
include("martensen/system_geometry.jl")
include("martensen/singularities.jl")
include("martensen/system_matrices.jl")
include("martensen/solve.jl")
include("martensen/post_process.jl")

# Hess-Smith (Educational)
include("hess_smith/method.jl")
include("hess_smith/panel_geometry.jl")
# include("hess_smith/system_geometry.jl")
# include("hess_smith/singularities.jl")
# include("hess_smith/system_matrices.jl")
# include("hess_smith/solve.jl")
# include("hess_smith/post_process.jl")

##### ----- CORE FUNCTIONALITY ----- #####

# Convenience Functions
include("convenience_functions.jl")

# Dispatch Functions
include("universal_dispatch.jl")

# Utility Functions
include("universal_utilities.jl")
include("universal_geometry_utilities.jl")

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

##### ----- TYPES ----- #####

export Mfoil, Xfoil, Lewis, Martensen#, HessSmith

##### ----- FUNCTIONS ----- #####

# Convenience Functions
export analyze

end
