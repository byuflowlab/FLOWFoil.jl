@testset "4 Panel Test Open TE" begin
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.5; 0.0; 0.5; 0.01]
    mesh = FLOWFoil.generate_mesh(x, y)
    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = 10.0
    re = rho * Vinf * chord / mu
    mach = 0.0
    alpha = 0.0

    coeffs = FLOWFoil.assemble_vortex_coefficients(mesh)

    @test !isapprox(LinearAlgebra.det(coeffs), 0.0)
end
