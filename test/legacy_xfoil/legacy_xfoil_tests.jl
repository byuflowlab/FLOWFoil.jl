@testset "LegacyXfoil Tests" begin

    # extract geometry
    x = Float64[]
    y = Float64[]

    f = open("test\\legacy_xfoil\\naca2412.dat", "r")

    for line in eachline(f)
        entries = split(chomp(line))
        push!(x, parse(Float64, entries[1]))
        push!(y, parse(Float64, entries[2]))
    end

    close(f)

    coordinates = hcat(x, y)

    # set operating conditions
    alpha = -9:1:13
    re = 1e5
    mach = 0.0
    n = length(alpha)

    method = LegacyXfoil(re, mach)
    
    # analyze with LegacyXfoil
    outputs = FLOWFoil.analyze_lxf(coordinates, alpha; method=method)

    # Using Xfoil example. See Xfoil documentation `Automated Angle of Attack Sweep`.
    clxf = [-0.650846,-0.662429,-0.635637,-0.565933,-0.479794,-0.390118,-0.301726,-0.150388,0.071766,	0.263166,	0.404422,	0.512138,0.611810,0.708466,0.804453,0.893460,0.969537,1.020497,1.087529,1.173720,1.242721,1.256900,1.197794]
    cdxf = [0.064745, 0.046357, 0.03411, 0.026936, 0.022554, 0.019767, 0.018318, 0.017026, 0.017390, 0.016828, 0.015968, 0.015591, 0.015749,0.016397, 0.017169,	0.018133,0.019458,	0.023512,	0.029120,0.035266,0.043276,0.054220,0.066990]
    cmxf = [-0.049415,-0.043245,-0.035047,-0.029727,-0.025883,-0.023503,-0.021741,-0.030251,-0.049186,-0.065553,-0.069629,-0.066432,-0.061841,-0.057255,-0.052536,-0.046931,-0.039801,-0.030233,-0.023964,-0.020845,-0.016251,-0.007123,0.003651]
    convxf = ones(Bool, n)

    @test isapprox(outputs.cl, clxf, atol=1e-5)
    @test isapprox(outputs.cd, cdxf, atol=1e-5)
    @test isapprox(outputs.cm, cmxf, atol=1e-5)
    @test outputs.converged == convxf
end