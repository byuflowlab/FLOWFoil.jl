
@testset "Post Processing Tests" begin

    # - Very Basic Test - #

    x, z = AirfoilTools.naca4()
    coordinates = [x z .+ 1.0]
    method = Lewis(; body_of_revolution=[false])

    # Generate Problem Object
    outputs = analyze(coordinates; method=method)

    #TODO: Need to figure out how to test this...
end
