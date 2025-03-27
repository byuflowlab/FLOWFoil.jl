@testset "Solve Tests" begin

    #---------------------------------#
    #             PLANAR              #
    #---------------------------------#
    # - Very Basic Test - #

    A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    b = [1.0; 2.0; 3.0]

    # Solve Linear System
    solution = solve(InviscidSystem(A, b, [3]))

    @test all(solution.x .== [1.0; 2.0; 3.0])
end
