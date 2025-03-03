using FLOWFoil
using Test
using LinearAlgebra
using Xfoil
using ForwardDiff
using ReverseDiff

# include("martensen/panel_geometry_tests.jl")
# include("martensen/system_geometry_tests.jl")
# include("martensen/system_matrix_tests.jl")

dir = "martensen/"
for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
    include(joinpath(dir, file))
end

# TODO: update lewis tests
# dir = "lewis/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end

# dir = "hess_smith/"
# for file in filter(f -> endswith(f, "tests.jl"), readdir(dir))
#     include(joinpath(dir, file))
# end
#
# # Mfoil (Xfoil)
# include.(filter(contains(r".jl$"), readdir("mfoil"; join=true)))
