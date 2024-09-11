"""
    solve_inviscid(inviscid_system, mesh)

Solve the inviscid_system for the vortex and streamfunction strengths.

Outputs the InviscidSolution object which contains the inviscid_system.

**Arguments:**
- `inviscid_system::InviscidSystem` : inviscid_system to solve.
- `mesh::Mesh` : Mesh defining geometry (to put into solution object)

**Returns:**
 - `solution::InviscidSolution`

"""
function solve_inviscid(inviscid_system)

    # Solve System
    if ndims(inviscid_system.b) > 1

        # initialize output
        x = similar(inviscid_system.b)

        # Pre-compute factorization
        Afact = LinearAlgebra.factorize(ImplicitAD.fd_value(inviscid_system.A))

        # loop through dimensions, reusing the factorization of A
        for i in 1:size(inviscid_system.b)[2]
            x[:, i] = ImplicitAD.implicit_linear(
                inviscid_system.A, inviscid_system.b[:, i]; Af=Afact
            )
        end

    else
        # Just solve like normal
        x = ImplicitAD.implicit_linear(inviscid_system.A, inviscid_system.b)
        # x = inviscid_system.A \ inviscid_system.b
    end

    return InviscidSolution(x, inviscid_system)
end
