function solve(::Lewis, system_matrices)

    # # Solve System
    if ndims(system_matrices.b) > 1
        # initialize output
        x = similar(system_matrices.b)
        
        # Pre-compute factorization
        Afact = LinearAlgebra.factorize(ImplicitAD.fd_value(system_matrices.A))
        # loop through dimensions, reusing the factorization of A
        for i in 1:size(system_matrices.b, 2)
            x[:, i] = ImplicitAD.implicit_linear(
            system_matrices.A, system_matrices.b[:, i]; Af=Afact
            )
        end
    else
        # Just solve like normal
        x = ImplicitAD.implicit_linear(system_matrices.A, system_matrices.b)
    end

    return x
end
