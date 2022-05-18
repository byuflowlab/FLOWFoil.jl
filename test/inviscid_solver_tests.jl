
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

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem)
    gamma = amat \ psi_inf
    gammatot = [gamma[i, 1] * cosd(alpha) + gamma[i, 2] * sind(alpha) for i in 1:length(x)]

    #TODO: figure out how to test this.
end

@testset "4 Panel Test Closed TE" begin
    x = [1.0; 0.5; 0.0; 0.5; 1.0]
    y = [0.0; -0.5; 0.0; 0.5; 0.0]
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
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem)
    gamma = amat \ psi_inf

    #TODO: figure out how to test this.
end

@testset "NACA Airfoil" begin
    x, yu, yl = FLOWFoil.naca4()
    x = [x[end:-1:2]; x[1:end]]
    y = [yl[end:-1:2]; yu[1:end]]
    mesh1 = FLOWFoil.generatemesh(x, y)
    meshsystem = FLOWFoil.MeshSystem([mesh1], [1.0], [0.0], [[0.0; 0.0]])
    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = 10.0
    re = rho * Vinf * chord / mu
    mach = 0.0
    alpha = 0.0
    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    coeffs = FLOWFoil.assemblevortexcoefficients(meshsystem)

    # Check to make sure coefficient matrix is not singular
    #    @test !isapprox(LinearAlgebra.det(coeffs), 0.0)

    #   @testset "Joukowsky Airfoil Velocity and Pressure Distributions" begin
    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem)
    gamma = amat \ psi_inf
    #  end
end

@testset "Joukowsky Airfoil" begin
    center = [-0.1; 0.1]
    R = 0.5
    alpha = 0.0
    N = 360
    U = 10.0

    x, y, vmag, cpdist = FLOWFoil.joukowskysurface(center, R, alpha, U; N=N)

    mesh1 = FLOWFoil.generatemesh(x, y)
    meshsystem = FLOWFoil.MeshSystem([mesh1], [1.0], [0.0], [[0.0; 0.0]])

    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = U
    re = rho * Vinf * chord / mu
    mach = 0.0

    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    coeffs = FLOWFoil.assemblevortexcoefficients(meshsystem)

    # Check to make sure coefficient matrix is not singular
    #    @test !isapprox(LinearAlgebra.det(coeffs), 0.0)

    #   @testset "Joukowsky Airfoil Velocity and Pressure Distributions" begin
    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem)
    gamma = amat \ psi_inf
    #  end
end
