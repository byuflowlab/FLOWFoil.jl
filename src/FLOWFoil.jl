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
abstract type Outputs end
abstract type AuxOutputs end

# Mfoil (Xfoil)
include("mfoil/method.jl")

# Lewis (Axisymmetric)
include("lewis/method.jl")
include("lewis/geometry_utils.jl")
include("lewis/panel_geometry.jl")
include("lewis/system_geometry.jl")
include("lewis/singularities.jl")
include("lewis/system_matrices.jl")
include("lewis/solve.jl")
include("lewis/post_process.jl")

# Martensen (Periodic)

# Hess-Smith (Educational)

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

export Mfoil, Xfoil, Lewis#, Martensen

##### ----- FUNCTIONS ----- #####

# Convenience Functions
export analyze

end
