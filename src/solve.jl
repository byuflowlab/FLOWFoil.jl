#=
Solve Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                          INVISCID SOLVER                           #
#                                                                    #
######################################################################

"""
    InviscidSolution{TM,TF,TD}

**Fields:**
 - `x::Array{Float}` : solution value(s) at each airfoil node.
 - `mesh::Mesh` : Mesh object describing airfoil nodes etc.
 - `system::InviscidSystem` : system object.
"""
struct InviscidSolution{TF,TS<:System}
    x::TF
    system::TS
end

