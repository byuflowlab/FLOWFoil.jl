using FLOWFoil
import FLOWFoil.AirfoilTools as at

using Test

using LinearAlgebra

using ForwardDiff
using FiniteDiff
using ImplicitAD

import NeuralFoil as nf
import Xfoil as xf

@info("Inviscid Mfoil Tests")
dir = joinpath(@__DIR__, "mfoil")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Axisymmetric Panel Method Tests")
dir = joinpath(@__DIR__, "lewis")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Periodic Panel Method Tests")
dir = joinpath(@__DIR__, "martensen")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("Hess Smith Tests")
dir = joinpath(@__DIR__, "hess_smith")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("LegacyXfoil Tests")
dir = joinpath(@__DIR__, "legacy_xfoil")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

@info("NeuralFoil Tests")
dir = joinpath(@__DIR__, "neural_foil")
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

