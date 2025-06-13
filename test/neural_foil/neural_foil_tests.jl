@testset "NeuralFoil Tests" begin

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
    flow_angles = -9:1:13
    reynolds = 1e6
    mach = 0.0

    nf_outputs = nf.get_aero_from_coordinates(coordinates, flow_angles, reynolds; model_size="xlarge")
    
    outputs = FLOWFoil.analyze_nf([x y], flow_angles; method=NeuralFoil())

    clnf = nf_outputs.cl
    cdnf = nf_outputs.cd
    cmnf = nf_outputs.cm
    acnf = nf_outputs.analysis_confidence

    @test isapprox(outputs.cl, clnf[:, 1], atol=1e-5)
    @test isapprox(outputs.cd, cdnf[:, 1], atol=1e-5)
    @test isapprox(outputs.cm, cmnf[:, 1], atol=1e-5)
    @test isapprox(outputs.analysis_confidence, acnf[:, 1], atol=1e-5)
end