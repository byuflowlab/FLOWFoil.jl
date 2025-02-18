function solve(::Martensen, system_matrices)
    #solve system
    gamma = ImplicitAD.linear_solve(system_matrices.A, system_matrices.b)

    #add on last strength based on kutta condition method
    return [gamma; -gamma[1, :]']
end
