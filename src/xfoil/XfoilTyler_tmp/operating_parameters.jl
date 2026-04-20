abstract type InviscidProblem end

struct HessSmith <: InviscidProblem end
struct XFOIL <: InviscidProblem end

struct OperatingParameters{TF}
    numpanels_airfoil::Int64
    blunt_TE::Bool
    Re::TF
    Mach::TF
    rho::TF
    mu::TF
    Vinf::TF
    chord::TF
    alpha::AbstractVector{TF} # in radians
    inviscidproblem::InviscidProblem
    
    # Viscous-specific
    viscous::Bool
    numpanels_wake::Int64
    ncrit::Float64
    wake_length::Float64 #percentage of chord to extend the wake
    etol::Float64
    maxiters::Int64
    γair::Float64
    GA::Float64
    GB::Float64
    GC::Float64
    H0::Float64         #stagnation enthalpy
    mu0::Float64        #stagnation viscosity
    rho0::Float64       #stagnation density
    KTb::Float64        #Karman-Tsien beta
    KTl::Float64        #Karman-Tsien lambda
    Cτ0::Float64
    Eτ0::Float64
    wake_gap::AbstractVector # wake TE gap
    use_native_solver_local::Bool
    transition_method::Int64 #1 is normal like mfoil (sum(lam, turb) in global solve), 2 is nested laminar solve, 3 is laminar & turbulent residuals in global solve #TODO: redo these numbers - may just need 1 and 2
    xft_l::TF
    xft_u::TF
    sft_l::Vector{TF}
    sft_u::Vector{TF}
    forced_l::Bool
    forced_u::Bool
    verbose::Bool
    
    # don't have defaults or not user-defined
    transition_lower::Vector{TF}
    transition_upper::Vector{TF}
    sinalpha::Vector{TF}
    cosalpha::Vector{TF} 
    separation_lower::Vector{TF}
    separation_upper::Vector{TF}
end
Base.eltype(::OperatingParameters{TF}) where TF = TF

function OperatingParameters(;numpanels_airfoil=200, blunt_TE=true,
                        Re=1e6, Mach=0.0, rho=1.0, mu=1e-6, Vinf=1.0, chord=1.0, 
                        alpha=nothing, inviscidproblem=XFOIL(), viscous=true, 
                        numpanels_wake=Int(ceil((numpanels_airfoil+1)/10 + 10) - 1),
                        ncrit=9.0, wake_length=1.0, etol=1e-10, maxiters=100, γair=1.4, GA=6.7, GB=0.75, GC=18.0, 
                        H0=0.0, mu0=0.0, rho0=1.0, KTb = 0.0, KTl = 0.0,
                        Cτ0=1.8, Eτ0=3.3, wake_gap=zeros(numpanels_wake+1), 
                        use_native_solver_local=true, transition_method=1, 
                        xft_l=chord, xft_u=chord, verbose=false, cosalpha=nothing, sinalpha=nothing) 

                        @assert (!isnothing(alpha) || !isnothing(cosalpha) || !isnothing(sinalpha)) "alpha or cosalpha or sinalpha must be given"

                        # if the user provides alpha, use that to get cosalpha and sinalpha, even if those are also provided
                        if !isnothing(alpha)
                            alpha = typeof(alpha)<:AbstractVector ? alpha : [alpha]
                            sinalpha = sin.(alpha)
                            cosalpha = cos.(alpha)
                        else #alpha was not specified, so use atan to extract it (though I shouldn't actually need it) #TODO: remove alpha?
                            if isnothing(cosalpha) #only sinalpha is given
                                sinalpha = typeof(sinalpha)<:AbstractVector ? sinalpha : [sinalpha]
                                alpha = asin(sinalpha)
                                cosalpha = cos.(alpha)
                            elseif isnothing(sinalpha) #only cosalpha is given
                                cosalpha = typeof(cosalpha)<:AbstractVector ? cosalpha : [cosalpha]
                                alpha = acos(cosalpha)
                                sinalpha = sin.(alpha)
                            else #both are given, this is ideal
                                sinalpha = typeof(sinalpha)<:AbstractVector ? sinalpha : [sinalpha]
                                cosalpha = typeof(cosalpha)<:AbstractVector ? cosalpha : [cosalpha]
                                alpha = atan.(sinalpha, cosalpha)
                            end
                        end

                        TF = promote_type(eltype(alpha), eltype(cosalpha), eltype(sinalpha), 
                                            typeof(chord))
                        alpha = convert(Vector{TF}, alpha)
                        cosalpha = convert(Vector{TF}, cosalpha)
                        sinalpha = convert(Vector{TF}, sinalpha)

                        Re     = convert(TF, Re)
                        Mach   = convert(TF, Mach)
                        rho    = convert(TF, rho)
                        mu     = convert(TF, mu)
                        Vinf   = convert(TF, Vinf)
                        chord  = convert(TF, chord)
                        xft_l = convert(TF, xft_l)
                        xft_u = convert(TF, xft_u)
                        sft_l = ones(TF, 1)
                        sft_u = ones(TF, 1)
                        forced_l = xft_l < chord
                        forced_u = xft_u < chord
                        transition_lower = ones(TF, 1)
                        transition_upper = ones(TF, 1)
                        separation_lower = ones(TF, 1)
                        separation_upper = ones(TF, 1)
                        
                        return OperatingParameters(
                            numpanels_airfoil, blunt_TE, Re, Mach, rho, mu, Vinf, chord, 
                            alpha,
                            inviscidproblem, viscous, numpanels_wake, ncrit, wake_length, 
                            etol, maxiters, γair, GA, GB, GC, 
                            H0, mu0, rho0, KTb, KTl, Cτ0, Eτ0, 
                            wake_gap, use_native_solver_local, transition_method, 
                            xft_l, xft_u, sft_l, sft_u, forced_l, forced_u,
                            verbose,
                            transition_lower, transition_upper,
                            sinalpha,
                            cosalpha,
                            separation_lower, separation_upper
                        )
end