
@testset "4 Panel Test Open TE" begin
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [-0.01; -0.25; 0.0; 0.25; 0.01]
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

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem, freestream)
    gamma = amat \ psi_inf
end

@testset "4 Panel Test Closed TE" begin
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [0.0; -0.25; 0.0; 0.25; 0.0]
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

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem, freestream)
    gamma = amat \ psi_inf
end

@testset "Joukowsky Airfoil" begin
    center = [-0.1; 0.1]
    R = 1.0
    alpha = 0.0
    U = 10.0

    x, y, vmag, cp = FLOWFoil.joukowskysurface(center, R, alpha, U; N=100)

    x, yu, yl = naca4(2.0, 4.0, 12.0, 80; bluntTE=false)
    x = [x[end:-1:1]; x[2:end]]
    y = [yl[end:-1:1]; yu[2:end]]
    mesh1 = FLOWFoil.generatemesh(x, y)
    meshsystem = FLOWFoil.MeshSystem([mesh1], [1.0], [0.0], [[0.0; 0.0]])
    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = 10.0
    re = rho * Vinf * chord / mu
    mach = 0.0
    alpha = 5.0
    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    coeffs = FLOWFoil.assemblevortexcoefficients(meshsystem)

    @test !isapprox(LinearAlgebra.det(coeffs), 0.0)

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem, freestream)
    gamma = amat \ psi_inf
end
