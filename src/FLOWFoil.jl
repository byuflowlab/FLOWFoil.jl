module FLOWFoil

#---------------------------------#
#           DEPENDENCIES          #
#---------------------------------#
using LinearAlgebra
using FLOWMath
using ImplicitAD
using SpecialFunctions
import NeuralFoil as nf
import Xfoil as xf

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

# NeuralFoil Translation
include("neural_foil/method.jl")

# LegacyXfoil
include("legacy_xfoil/method.jl")

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

export Mfoil, Xfoil, Lewis, Martensen, LegacyXfoil, NeuralFoil
export InviscidOutputs, LegacyXFOutputs#, NeuralOutputs

##### ----- FUNCTIONS ----- #####

# Convenience Functions
export analyze

end
