function solve(::Martensen, system_matrices)
    #solve system
    TF = eltype(system_matrices.b)
    n = length(system_matrices.b[:, 1])
    x = Array{TF, 2}(undef, n, 2) #x is the solution to the system of linear equations
    
    x[:, 1] = system_matrices.A \ system_matrices.b[:, 1]
    x[:, 2] = system_matrices.A \ system_matrices.b[:, 2]

    return x
end

######################################################################
#                                                                    #
#                        FROM CHATGPT TRANSLATION OF LEWIS                         #
#                                                                    #
######################################################################

function solution()
    # Declare variables
    i = 0
    m = 0  # Assuming m is defined somewhere else or provided
    gamma1 = 0.0
    gamma2 = 0.0
    pitch = 0.0
    beta1 = 0.0
    beta2 = 0.0
    betainf = 0.0
    W1 = 0.0
    Winf = 0.0
    Uinf = 0.0
    Vinf = 0.0
    ans1 = zeros(Float64, m)  # Assuming ans1 is defined elsewhere
    ans2 = zeros(Float64, m)  # Assuming ans2 is defined elsewhere
    Cp = 0.0
    gamma = 0.0
    chord = 0.0

    # Calculate k1, k2, beta2, and betainf
    k1 = (1.0 - gamma2 / 2.0 / pitch) / (1.0 + gamma2 / 2.0 / pitch)
    k2 = gamma1 / pitch / (1.0 + gamma2 / 2.0 / pitch)
    beta2 = atan(k1 * sin(beta1) / cos(beta1) - k2)
    betainf = atan(0.5 * (sin(beta1) / cos(beta1) + sin(beta2) / cos(beta2)))

    # Calculate Winf, Uinf, and Vinf
    Winf = W1 * cos(beta1) / cos(betainf)
    Uinf = Winf * cos(betainf)
    Vinf = Winf * sin(betainf)

    # Print elements and Cp
    println("Element ", rpad("Cp", 9))
    println()

    for i in 1:m
        ans = Uinf * ans1[i] + Vinf * ans2[i]
        Cp = 1.0 - sqrt(ans / W1)^2
        println(rpad(i, 4), rpad(Cp, 16, ":6"))
    end

    # Calculate gamma and other results
    gamma = Uinf * gamma1 + Vinf * gamma2
    println("Predicted circulation = ", lpad(gamma, 10, ":6"))
    println("Cl = ", lpad(2.0 * gamma / (Winf * chord), 10, ":6"))
    return println("Beta2 = ", lpad(beta2 * 180 / pi, 10, ":6"))
end

