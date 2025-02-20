function solve(::HessSmith, system_matrices)
    return ImplicitAD.linear_solve(system_matrices.A, system_matrices.b)
end
