# Start with conformal mapping solution
function conformal_map(U=10.0, rho=1.225; alpha=0.0)

    # Karman Trefftz
    x, y, th, R, beta = FLOWFoil.rawjoukowsky([-0.1 0.1])

    # analytic solutions
    gamma = -4.0 * pi * R * U * sin(beta + alpha)
    vth = -(2.0 * U * sin.(th) .+ gamma / (2 * pi * R))
    pth = rho * U^2 * cos.(2.0 .* th)

    z = R .* exp.(im .* collect(th))
    wz = U .* (z .+ R^2.0 ./ z) .+ im .* (gamma .* log.(z) ./ (2.0 .* pi))
    dwz = U .* (1.0 .- R .^ 2.0 ./ z .^ 2.0) + im .* (gamma ./ (2.0 .* pi .* z))
    vx = real.(dwz)
    vy = imag.(dwz)

    # shift everything around so x,y coordinates are usable for panel method
    # find maximum x
    x = reverse(x)
    y = reverse(y)
    th = reverse(th)
    vth = reverse(vth)
    pth = reverse(pth)

    return vth, pth, x, y, th
end

@testset "4 Panel Test" begin
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
    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem, freestream)
end

@testset "Joukowsky Airfoil" begin
    vth, pth, x, y, th = conformal_map(; alpha=5.0 * pi / 180.0)
    mesh = FLOWFoil.generatemesh(x, y)
    meshsystem = FLOWFoil.MeshSystem([mesh], [1.0], [0.0], [[0.0; 0.0]])
    chord = maximum(x) - minimum(x)
    rho = 1.225
    mu = 1.8e-5
    Vinf = 10.0
    re = rho * Vinf * chord / mu
    mach = 0.0
    freestream = FLOWFoil.Freestream([re], [rho], [mu], [mach], [alpha])

    amat = FLOWFoil.assemblematrixa(meshsystem)
    psi_inf = FLOWFoil.assembleboundaryconditions(meshsystem, freestream)
end
