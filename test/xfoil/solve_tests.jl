@testset "Solve Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    b = [1.0; 2.0; 3.0]

    # Solve Linear System
    solution = FLOWFoil.solve_inviscid(Xfoil(), (;A, b))

    @test all(solution .== [1.0; 2.0; 3.0])
end
