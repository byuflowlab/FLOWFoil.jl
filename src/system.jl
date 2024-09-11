#=

Inviscid System Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                              GENERAL                               #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

abstract type System end

"""
    InviscidSystem

**Fields:**
 - `A::Array{Float,2}` : Coefficient Matrix on Left Hand Side.
 - `b::Array{Float,2}` : Boundary Condition Coefficient Vector on Right Hand Side.
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
struct InviscidSystem{TA,TB,TI} <: System
    A::TA
    b::TB
    Ns::TI
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

######################################################################
#                                                                    #
#                               PLANAR                               #
#                                                                    #
######################################################################

#= NOTE:
The system assembly here is based on the Xfoil implementation, which does not strictly fit in the overall structure (it solves based on stream functions rather than potentials).
Likely, this will be moved elsewhere as other methods are developed.
=#
######################################################################
#                                                                    #
#                            AXISYMMETRIC                            #
#                                                                    #
######################################################################

######################################################################
#                                                                    #
#                             PERIODIC                               #
#                                                                    #
######################################################################

