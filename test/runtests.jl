using FLOWFoil
using Test
using LinearAlgebra
using Xfoil
using ForwardDiff
using ReverseDiff

include("problem_tests.jl")
# include("autodiff_tests.jl")

# Mfoil (Xfoil)
include.(filter(contains(r".jl$"), readdir("mfoil"; join=true)))

# Lewis (Axisymmetric)
include.(filter(contains(r".jl$"), readdir("lewis"; join=true)))

# Martensen (Periodic)
include.(filter(contains(r".jl$"), readdir("martensen"; join=true)))

# Hess-Smith (Educational)
include.(filter(contains(r".jl$"), readdir("hess_smith"; join=true)))
