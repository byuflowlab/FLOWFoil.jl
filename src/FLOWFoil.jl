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

##### ----- CORE FUNCTIONALITY ----- #####

# Convenience Functions
include("convenience_functions.jl")

# Dispatch Functions
include("universal_dispatch.jl")

# Utility Functions
include("universal_utilities.jl")

##### ----- Methods ----- #####
# Mfoil (Xfoil)
include.(filter(contains(r".jl$"), readdir("mfoil"; join=true)))

# Lewis (Axisymmetric)
include.(filter(contains(r".jl$"), readdir("lewis"; join=true)))

# Martensen (Periodic)
include.(filter(contains(r".jl$"), readdir("martensen"; join=true)))

# Hess-Smith (Educational)
include.(filter(contains(r".jl$"), readdir("hess_smith"; join=true)))

#---------------------------------#
#             EXPORTS             #
#---------------------------------#

##### ----- TYPES ----- #####

export Mfoil, Xfoil, Lewis#, Martensen

##### ----- FUNCTIONS ----- #####

# Convenience Functions
export solve
