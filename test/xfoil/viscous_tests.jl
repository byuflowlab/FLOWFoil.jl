#=
Tests for the coupled viscous + inviscid Xfoil-style solve.

These exercise the full pipeline that runs when `Xfoil(viscous=true)` is
passed to `analyze`. Closures and station residuals are tested implicitly
via the global Newton in the convergence tests (the global solver invokes
them thousands of times per case).

NACA 4-digit airfoils are constructed with `blunt_trailing_edge=true` —
mfoil and the original Drela formulation are designed around blunt TE
airfoils; sharp-TE viscous is out of scope.
=#

@testset "Viscous: NACA 2412 @ α=5 deg (canonical)" begin
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [5.0];
        method=Xfoil(
            viscous=true,
            Re=1e6,
            Mach=0.0,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            verbose=false,
        ),
    )

    @test outputs.converged[1]
    # Reference values: mfoil-style coupled solve on the
    # same blunt-TE NACA 2412 give cl ≈ 0.81, cd ≈ 0.008, transition near
    # 30-35% chord on the upper surface and no transition on the lower.
    @test 0.78 < outputs.cl[1] < 0.84
    @test 0.006 < outputs.cd[1] < 0.010
    @test 0.25 < outputs.transition_upper[1] < 0.40
    @test outputs.transition_lower[1] ≈ 1.0
end

@testset "Viscous: NACA 2412 @ α=0 deg" begin
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [0.0];
        method=Xfoil(
            viscous=true,
            Re=1e6,
            Mach=0.0,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            verbose=false,
        ),
    )

    @test outputs.converged[1]
    # cambered NACA 2412 still has cl > 0 at α=0
    @test 0.18 < outputs.cl[1] < 0.30
    @test 0.004 < outputs.cd[1] < 0.010
end

@testset "Viscous: multi-alpha sweep" begin
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [0.0, 2.0, 5.0];
        method=Xfoil(
            viscous=true,
            Re=1e6,
            Mach=0.0,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            verbose=false,
        ),
    )

    @test all(outputs.converged)
    # cl rises monotonically with α
    @test outputs.cl[1] < outputs.cl[2] < outputs.cl[3]
    # Transition moves forward with α (for ncrit = 9, blunt-TE NACA 2412)
    @test outputs.transition_upper[1] > outputs.transition_upper[3]
end

@testset "Viscous: forced transition" begin
    # Force transition at 50% chord on both surfaces.
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [0.0];
        method=Xfoil(
            viscous=true,
            Re=1e5,
            Mach=0.0,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            xft_l=0.5,
            xft_u=0.5,
            verbose=false,
        ),
    )

    @test outputs.converged[1]
    # Forced transition at x/c=0.5: actual transition should be at or near 0.5
    @test outputs.transition_upper[1] ≤ 0.55
    @test outputs.transition_lower[1] ≤ 0.55
end

@testset "Viscous: compressible NACA 2412 @ M=0.3, α=5 deg" begin
    # Karman-Tsien compressibility path.
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [5.0];
        method=Xfoil(
            viscous=true,
            Re=1e6,
            Mach=0.3,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            verbose=false,
        ),
    )

    @test outputs.converged[1]
    @test isapprox(outputs.cl[1], 0.8686; atol=1e-3)
    @test isapprox(outputs.cd[1], 0.008266; atol=1e-4)
    @test isapprox(outputs.transition_upper[1], 0.3164; atol=1e-3)
    @test outputs.transition_lower[1] ≈ 1.0
end

@testset "Viscous: compressible NACA 2412 @ M=0.5, α=3 deg" begin
    # Higher-Mach compressible case to exercise the KTl/KTb branches further.
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [3.0];
        method=Xfoil(
            viscous=true,
            Re=1e6,
            Mach=0.5,
            ncrit=9.0,
            etol=1e-10,
            maxiters=50,
            verbose=false,
        ),
    )

    @test outputs.converged[1]
    @test isapprox(outputs.cl[1], 0.7114; atol=1e-3)
    @test isapprox(outputs.cd[1], 0.007118; atol=1e-4)
    @test isapprox(outputs.transition_upper[1], 0.4415; atol=1e-3)
    # Lower side stays essentially laminar (any transition is within the last ~1% chord)
    @test outputs.transition_lower[1] > 0.95
end

@testset "Viscous: output structure" begin
    x, y = at.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)
    coordinates = [x y]
    outputs = analyze(
        coordinates,
        [3.0];
        method=Xfoil(viscous=true, Re=1e6, etol=1e-10, maxiters=50, verbose=false),
    )

    # Inviscid intermediates are exposed alongside viscous results
    @test outputs.gamma_ref isa AbstractMatrix
    @test size(outputs.gamma_ref, 2) == 2          # cos α and sin α columns
    @test outputs.A isa AbstractMatrix
    @test size(outputs.A, 1) == size(outputs.A, 2)

    # Per-α scalars indexed by α
    @test length(outputs.cl) == 1
    @test length(outputs.cd) == 1
    @test length(outputs.cm) == 1
    @test length(outputs.converged) == 1
    @test length(outputs.transition_upper) == 1
    @test length(outputs.transition_lower) == 1

    # Per-α detailed BL state
    @test length(outputs.per_alpha) == 1
    pa = outputs.per_alpha[1]
    @test pa.deltas_u isa AbstractVector
    @test pa.thetas_l isa AbstractVector
    @test pa.state3s_w isa AbstractVector
end
