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
    AxisymmetricPolar{TF,TA}

**Fields:**
- `surface_velocity::Matrix{Float}` : surface velocity on each panel
- `surface_pressure::Matrix{Float}` : surface pressure coefficient on each panel
- `xsmooth::Matrix{Float}` : x-values associated with smoothed surface distributions
"""
struct AxisymmetricPolar{TF} <: Polar
    surface_velocity::Matrix{TF}
    surface_pressure::Matrix{TF}
    xsmooth::Matrix{TF}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#


######################################################################
#                                                                    #
#                     PERIODIC POST PROCESSING                       #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    PeriodicPolar{TF,TA}

**Fields:**
- `surface_velocity::Matrix{Float}` : surface velocity on each panel
- `surface_pressure::Matrix{Float}` : surface pressure coefficient on each panel
- `xsmooth::Matrix{Float}` : x-values associated with smoothed surface distributions
"""
struct PeriodicPolar{TF} <: Polar
    surface_velocity::Array{TF,3}
    surface_pressure::Array{TF,3}
    xsmooth::Array{TF,3}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#


