#=
Compose Type Definitions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                          VISCOUS TYPES                             #
#                                                                    #
######################################################################



# """
# **Fields:**
# """
# struct ViscousSystem{TPa,TPr,TF}
#     parameters::TPa
#     properties::TPr
#     vcoeffmat::Array{TF}
#     scoeffmat::Array{TF}
#     bccoeffvec::Array{TF}
# end


# """
#     WakeMesh{TF}

# **Fields:**
#  - `wake_nodes::Array{Float,2}` : x,y wake panel node locations.
#  - `wake_midpoints::Array{Float,2}` : x,y wake panel center point locations.
# """
# struct WakeMesh{TF}
#     wake_nodes::Array{TF}
#     wake_midpoints::Array{TF}
# end

# """
# ViscousSolution{}

# **Fields:**
#  - `panelgammas::Array{Float,2}` : \$\\gamma_0\$ and \$\\gamma_{90}\$ values at each airfoil node.
#  - `panelsources::Array{Float,2}` : source values at each airfoil node.
#  - `wakesources::Array{Float,2}` : source values at each wake node.
#  - `bodystrength::Array{Float}` : \$\\Psi_0\$ (constant stream function) 0 and 90 values.
# """
# struct ViscousSolution{TF,TD}
#     panelgammas::TF
#     panelsources::TF
#     wakesources::TF
#     bodystrength::TF
# end

#"""
#    Properties{TF}

#Thermodynamic properties for the viscous solution

#**Fields:**
# - `machinf::Float` : freestream mach number
# - `KTbeta::Float` : Karman-Tsien beta
# - `KTlambda::Float` : Karman-Tsien lambda
# - `H0::Float` : stagnation enthalpy
# - `sonic_cp::Float` : sonic cp
# - `rho0::Float` : stagnation density
# - `mu0::Float` : stagnation dynamic viscosity
#"""
#struct Properties{TF}
#    machinf::TF
#    KTbeta::TF
#    KTlambda::TF
#    H0::TF
#    sonic_cp::TF
#    rho0::TF
#    mu0::TF
#end

#"""
#    Parameters{TF}

#Solver Parameters.

#**Fields:**
# - `gamma_air::Float = 1.4` : ratio of specific heats for air
# - `eta_crit::Float = 9.0` : critical amplification factor
# - `eta_D::Float = 0.9` : wall/wake dissipation length ratio
# - `GA::Float = 6.7` : G - Beta locus A constant
# - `GB::Float = 0.75` : G - Beta locus B constant
# - `GC::Float = 18.0` : G - Beta locus C constant
# - `Klag::Float = 5.6` : shear lag constant
# - `Ctau::Float = 1.8` : shear stress initialization constant
# - `Etau::Float = 3.3` : shear stree initialization exponent
# - `rSu::Float = 0.35` : Sutherland temperature ratio
# - `fw::Float = 2.5` : wake gap continuation factor
# - `dw::Float = 1.0` : wake length, in airfoil chords
# - `epsilonw::Float = 1e-5` : first wake point offset, in airfoil chords
# - `iknowwhatimdoing::Bool` : boolean to silence warnings if you really know what you're doing.
# - `rhoinf::Float` : non-dimensional freestream density
# - `vinf::Float` : non-dimensional freestream velocity magnitude
# - `muinf::Float` : freestream dynamic viscosity
#"""
#struct Parameters{TF,TB}
#    gamma_air::TF
#    etacrit::TF
#    etaD::TF
#    GA::TF
#    GB::TF
#    GC::TF
#    Klag::TF
#    Ctau::TF
#    Etau::TF
#    rSu::TF
#    fw::TF
#    dw::TF
#    epsilonw::TF
#    iknowwhatimdoing::TB
#    rhoinf::TF
#    vinf::TF
#    muinf::TF
#end

#"""
#    defaultparameters(;kwargs)

#Initialized parameters to defaults, but allows selective user override through keyword arguments.

#**Keyword Arguments:**
# - `gamma_air::Float` : ratio of specific heats for air
# - `eta_crit::Float` : critical amplification factor
# - `eta_D::Float` : wall/wake dissipation length ratio
# - `GA::Float` : G - Beta locus A constant
# - `GB::Float` : G - Beta locus B constant
# - `GC::Float` : G - Beta locus C constant
# - `Klag::Float` : shear lag constant
# - `Ctau::Float` : shear stress initialization constant
# - `Etau::Float` : shear stree initialization exponent
# - `rSu::Float` : Sutherland temperature ratio
# - `fw::Float` : wake gap continuation factor
# - `dw::Float` : wake length, in airfoil chords
# - `epsilonw::Float` : first wake point offset, in airfoil chords
# - `iknowwhatimdoing::Bool` : boolean to silence warnings if you really know what you're doing.
# - `rhoinf::Float` : non-dimensional freestream density
# - `vinf::Float` : non-dimensional freestream velocity magnitude
# - `muinf::Float` : freestream dynamic viscosity
#"""
#function defaultparameters(;
#    gamma_air=1.4,
#    eta_crit=9.0,
#    eta_D=0.9,
#    GA=6.7,
#    GB=0.75,
#    GC=18.0,
#    Klag=5.6,
#    Ctau=1.8,
#    Etau=3.3,
#    rSu=0.35,
#    fw=2.5,
#    dw=1.0,
#    epsilonw=1e-5,
#    iknowwhatimdoing=false,
#    rhoinf=1.0,
#    vinf=1.0,
#    muinf=0.0,
#)

#    #check if vinf and rhoinf are not default
#    if !iknowwhatimdoing && (vinf != 1.0 || rhoinf != 1.0)
#        @warn(
#            "vinf and/or rhoinf have been changed from their defaults. Do you really want to do that??\n\nTo silence this warning in the future, set the field `iknowwhatimdoing=true`."
#        )
#    end

#    return Parameters(
#        gamma_air, eta_crit, eta_D, GA, GB, GC, Klag, Ctau, Etau, rSu, fw, dw, epsilonw
#    )
#end
