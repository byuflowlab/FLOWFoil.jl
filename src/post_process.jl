#=
General Post Processing Types and Functions

Authors: Judd Mehr,
=#

abstract type Polar end

######################################################################
#                                                                    #
#                       PLANAR POST PROCESSING                       #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    PlanarPolar{TF}

Also used for Periodic (cascade) post processing.

**Fields:**
 - `lift::Float` : Lift Coefficient.
 - `drag::Float` : Total Drag Coefficient.
 - `pdrag::Float` : Pressure Drag Coefficient.
 - `idrag::Float` : Induced Drag Coefficient.
 - `moment::Float` : Moment Coefficient.
 - `surfacevelocity::Vector{Float}` : surface velocity distribution
 - `surfacepressure::Vector{Float}` : surface pressure distribution
"""
struct PlanarPolar{TF,TS<:Vector{TF}} <: Polar
    lift::TF
    drag::TF
    pdrag::TF
    idrag::TF
    moment::TF
    surface_velocity::TS
    surface_pressure::TS
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

######################################################################
#                                                                    #
#                    AXISYMMETRIC POST PROCESSING                    #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    AxiSymPolar{TF,TA}

**Fields:**
- `thrust::Float` : Thrust (or drag) of body
- `surface_velocity::Array{Float}` : surface velocity on each panel
- `surface_pressure::Array{Float}` : surface pressure coefficient on each panel
"""
struct AxiSymPolar{TF,TA} <: Polar
    thrust::TF
    surface_velocity::TA
    surface_pressure::TA
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#
