function solve(method::Xfoil, system_matrices)
    if method.viscous
        error("Viscous method not yet implemented.")
    else
        return solve_inviscid(method, system_matrices)
    end
end

"""
    solve_inviscid(system_matrices, mesh)

solve the system_matrices for the vortex and streamfunction strengths.

Outputs the InviscidSolution object which contains the system_matrices.

# Arguments
- `system_matrices::InviscidSystem` : system_matrices to solve.
- `mesh::Mesh` : Mesh defining geometry (to put into solution object)

# Returns
 - `x::inviscid solution, consisting of node vortex strengths`
"""
function solve_inviscid(method::Xfoil, system_matrices)

    # solve System
    if ndims(system_matrices.b) > 1

        # initialize output
        x = similar(system_matrices.b)

        # Pre-compute factorization
        Afact = LinearAlgebra.factorize(ImplicitAD.fd_value(system_matrices.A))

        # loop through dimensions, reusing the factorization of A
        for i in 1:size(system_matrices.b)[2]
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
