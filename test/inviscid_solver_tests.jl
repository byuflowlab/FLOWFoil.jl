@testset "4 Panel Test Open TE" begin
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01]
    mesh = FLOWFoil.generatemesh(x, y)
    meshsystem = FLOWFoil.MeshSystem([mesh], [1.0], [0.0], [[0.0; 0.0]])
    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = 10.0
    re = rho * Vinf * chord / mu
    mach = 0.0
    alpha = 0.0
    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    coeffs = FLOWFoil.assemblevortexcoefficients(meshsystem)

    @test !isapprox(LinearAlgebra.det(coeffs), 0.0)
end
