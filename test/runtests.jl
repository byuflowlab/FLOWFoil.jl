using FLOWFoil
import FLOWFoil.AirfoilTools as at
using Test
using LinearAlgebra
using Xfoil
using ForwardDiff
using FiniteDiff
using ImplicitAD

@info("Hess Smith Tests")
dir = "hess_smith/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Axisymmetric Panel Method Tests")
dir = "lewis/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Periodic Panel Method Tests")
dir = "martensen/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Inviscid Mfoil Tests")
dir = "mfoil/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("NeuralFoil Translation Tests")
dir = "neural_foil/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

