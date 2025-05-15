using FLOWFoil
using Test
using LinearAlgebra
using Xfoil
using ForwardDiff
using ReverseDiff
using ImplicitAD

# dir = "hess_smith/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end
#
# dir = "lewis/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end
#
# dir = "martensen/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end
#
# dir = "mfoil/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end

dir = "neural_foil/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

