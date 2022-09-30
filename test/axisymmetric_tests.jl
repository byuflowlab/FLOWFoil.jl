@testset "Body of Revolution" begin
    include("data/bodyofrevolutioncoords.jl")

    mesh = [FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=true)]

    problem = FLOWFoil.Problem(mesh; axisymmetric=true, viscous=false)

    inviscid_solution = FLOWFoil.solve(problem)

end
