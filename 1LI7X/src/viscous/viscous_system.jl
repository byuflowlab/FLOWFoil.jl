#=
Boudary Layer Integral Method with Transition Models

Authors: Andrew Ning

Date Started: 2 May 2023

Change Log:
=#


# ---------- parameters ----------
struct Params{TF}
    γ::TF  # ratio of specific heats
    ncrit::TF  # critical amplification factor in e^n method
    GA::TF  # G-β locus with A constant
    GB::TF  # G-β locus with B constant
    GC::TF  # G-β locus with C constant
    ηd::TF  # wall or wake dissipation length ratio
    Klag::TF  # shear lag constant
    Ctau::TF  # shear stress initialization constant
    Etau::TF  # shear stress initialization exponent
    rsu::TF  # sutherland temperature ratio
    fw::TF  # wake gap continuation factor
    dw::TF  # wake length, in airfoil chords
    ew::TF  # first wake oint offset in airfoil chords
end

# default values
Params() = Params(1.4, 9, 6.7, 0.75, 18.0, 0.9, 5.6, 1.8, 3.3, 0.35, 2.5, 1, 1e-5)

# ---------- compressibility --------------

# freestream parameters used in Karman-Tsien correction
function kt_params(Minf)
    β = sqrt(1 - Minf^2)
    λ = Minf^2 / (1 + β)^2
    return β, λ
end

# local Karman-Tsien velocity correction
kt_velocity(u, λ, Vinf) = u*(1 - λ) / (1 - λ*(u/Vinf)^2)
kt_cp(cp, λ, β) = cp / (β + λ*(1 + β)*cp/2)

# local speed of sound
speedsound(u, ht, γ) = sqrt((γ-1)*(ht - u^2/2))

# local mach number
mach(u, ht, γ) = u / speedsound(u, ht, γ)

# freestream total enthalpy
function totalenthapy(Vinf, Minf, γ)
    γm = γ-1
    factor = 1 + γm/2*Minf^2
    return Vinf^2 / (γm * Minf^2) * factor
end

# local density
function density(Me, ρinf, Minf, γ)
    γm2 = (γ-1)/2
    ex = 1/(γ - 1)
    return ρinf * (1 + γm2*Minf^2)^ex / (1 + γm2*Me^2)^ex
end

# ratio of local temperature to total temperature
tempratio(u, ht) = 1 - 0.5*u^2/ht

# subfunction in Sutherland's law
su_func(Tratio, rsu) = Tratio^1.5*(1 + rsu)/(Tratio + rsu)

# freestream parameters for Sutherland's
function su_params(Vinf, muinf, ht, rsu)
    Tr = tempratio(Vinf, ht)
    mu0 = muinf * su_func(Tr, rsu)
    return mu0
end

# local viscosity
function viscosity(u, mu0, rsu)
    Tr = tempratio(u, ht)
    mu = mu0 * su_func(Tr, rsu)
    return mu
end

# combine all the above for convenience
function compressibility(ue, θ, params)
    (; ht, λ, γ, ρinf, Minf, mu0, rsu) = params

    u = kt_velocity(u, λ, Vinf)
    M = mach(u, ht, γ)
    rho = density(M, ρinf, Minf, γ)
    mu = viscosity(u, mu0, rsu)
    Reθ = rho*u*θ/mu

    return u, M, Reθ
end

# set constants based on freestream inputs
function set_freestream(Minf, Reinf, Vinf, rhoinf, chord, γ)
    β, λ = kt_params(Minf)
    ht = totalenthapy(Vinf, Minf, γ)
    muinf = rhoinf * Vinf * chord / Reinf
    mu0 = su_params(Vinf, muinf, ht, rsu)

    return (β=β, λ=λ, ht=ht, mu0=mu0)
end


# -----------------------------------------

# ------ shape parameters --------

# all the shape factors
function shape(H, M, Reθ, bltype)

    # H_k
    Hk = (H - 0.29*M^2)/(1 + 0.113*M^2)

    if bltype == "wake"
        Hk = max(Hk, 1.00005)  # TODO
    else
        Hk = max(Hk, 1.05)
    end

    #H^{**}
    Hss = 0.064*M^2 / (Hk - 0.8)

    #H^*
    if bltype == "laminar"

        Hkt = Hk - 4.35
        if Hk < 4.35  # TODO: make these if/else statements smooth
            Hsl = (0.0111*Hkt^2 - 0.0278^Hkt^3) / (Hk - 1) + 1.528 - 0.002*(Hkt*Hk)^2
        else
            Hsl = 0.015*Hkt^2/Hk + 1.528
        end

        Hs = Hsl

    else  # turbulent

        H0 = min(3.0 + 400/Reθ, 4.0)  # TODO: replace with smooth min
        Hr = (H0 - Hk) / (H0 - 1)
        Reθt = max(Reθ, 200)  # TODO: replace with smooth max
        AH = Hk - H0 + 4/log(Reθt)
        if Hk < H0
            Hsti = 1.5 + 4/Reθt + 1.5*(0.5 - 4/Reθt)*Hr^2/(Hk + 0.5)
        else
            Hsti = 1.5 + 4/Reθt + (Hk - H0)^2*(0.007*log(Reθt)/AH^2 + 0.015/Hk)
        end
        Hst = (Hsti + 0.028*M^2) / (1 + 0.014*M^2)

        Hs = Hst
    end

    return Hk, Hss, Hs
end




# ------- skin friction  ----------------

# laminar skin friction
function cf_lam(Reθ, Hk)

    if Hk < 5.5  # TODO: make these if/else statements smooth
        num = 0.0727*(5.5 - Hk)^3 / (Hk + 1) - 0.07
    else
        num = 0.015*(1 - 1/(Hk - 4.5))^2 - 0.07
    end

    return num/Reθ
end

# turbulent skin friction
function cf_turb(Reθ, Hk, M, γ)

    Acf = -1.33*Hk
    if Acf < -17  # TODO
        Acf = -20.0 + 3*exp((Acf + 17)/3.0)
    end

    Fc = sqrt(1 + 0.5*(γ - 1)*M^2)
    Bcf = max(log10(Reθ/Fc), 1.303)  # TODO

    return 0.3*exp(Acf)*Bcf^(-1.74 - 0.31*Hk) + 0.00011*(tanh(4 - Hk/0.875) - 1)
end

# combine for convenience
function skinfriction(Reθ, Hk, M, γ, bltype)

    if bltype == "laminar"
        return cf_lam(Reθ, Hk)

    elseif bltype == "turbulent"
        return cf_turb(Reθ, Hk, M, γ)

    else  # wake
        return 0.0
    end

end

# ------- dissipation  ----------------

# laminar dissipation
function cdi_lam(Reθ, Hk)

    if Hk < 4  # TODO
        num = 0.00205*(4 - Hk)^5.5 + 0.207
    else
        num = -0.0016*(Hk - 4)^2 / (1 + 0.02*(Hk - 4)^2) + 0.207
    end

    return num/Reθ
end

# slip velocity
function slip_velocity(H, Hk, Hs, Gb, bltype)

    if bltype == "laminar"
        return 1.0  # not used
    end

    Us = 0.5*Hs*(1 - (Hk-1)/(H*Gb))

    if bltype == "wake"
        Us = min(Us, 0.99995)  # TODO
    else
        Us = min(Us, 0.98)
    end

    return Us
end

# sub function used in turbulent and wake coefficients
function cdi_sub(cτ, Reθ, Hs, Us)

    cdio = cτ * (0.995 - Us) * 2 / Hs  # outer
    cdis = 0.3 * (0.995 - Us)^2 / (Hs * Reθ)  # stress

    return cdio + cdis
end

# turbulent dissipation
function cdi_turb(cτ, Reθ, Hk, Hs, cf, Us)

    cdi = cdi_sub(cτ, Reθ, Hs, Us)

    # wall
    cdiw = 0.5 * cf * Us / Hs * (1 + tanh((Hk-1)*log(Reθ)/2.1))

    cdi += cdiw

    return min(cdi, cdi_lam(Reθ, Hk))  # TODO
end

# wake dissipation
function cdi_wake(cτ, Reθ, Hk, Hs, Us)

    cdi = dissipation_sub(cτ, Reθ, Hs, Us)

    cdilw = 2.2*(1 - 1/Hk)^2 / (Hk*Hs*Reθ)

    return min(cdi, cdilw)  # TODO

end

# combine for convenience
function dissipation(cτ, Reθ, H, Hk, Hs, cf, Gb, bltype)

    Us = slip_velocity(H, Hk, Hs, Gb, bltype)

    if bltype == "laminar"
        cdi = cdi_lam(Reθ, Hk)

    elseif bltype == "turbulent"
        cdi = cdi_turb(cτ, Reθ, Hk, Hs, cf, Us)

    else  # wake
        cdi = cdi_wake(cτ, Reθ, Hk, Hs, Us)
    end

    return Us, cdi
end
# ------------------------------------------

# ------ amplification: e^n method-------

# rate of increase in amplification factor
function amplification(θ, n, Reθ, Hk, params)

    epsn = 0.001*(1 + tanh(5*(n - params.ncrit)))
    Hhat = 1/(Hk - 1)
    L0 = 2.492*Hhat^0.43 + 0.7*(1 + tanh(14*Hhat - 9.24))
    sn = (log10(Reθ) - (L0 - 0.1)) / 0.2
    if sn < 0  # TODO
        rn = 0
    elseif sn > 1
        rn = 1
    else
        rn = 3*sn^2 - 2*sn^3
    end
    gn = 0.028/Hhat - 0.0345*exp(-(3.87*Hhat - 2.52)^2)
    fn = -0.05 + 2.7*Hhat - 5.5*Hhat^2 + 3*Hhat^3 + 0.1*exp(-20*Hhat)

    return (rn*fn*gn + epsn)/θ
end

# -----------------------------------------

# ------ lag -----------
# H_{kc}
function get_Hkc(Reθ, Hk, params, bltype)
    Gc = bltype == "wake" ? 0.0 : params.Gc
    Hkc = Hk - 1 - Gc/Reθ
    return Hkc
end

#u_q
function get_uq(δs, Reθ, Hk, cf, params, bltype)
    (; Ga, Gb, ηd) = params
    Hkc = get_Hkc(Reθ, Hk, params, bltype)
    uq = (0.5*cf - (Hkc/(Ga*ηd*Hk))^2) / (Gb * δs)
    return uq
end

# c_{\tau, eq}
function get_cτeq(Reθ, H, Hk, Hs, Us)
    (; Ga, Gb) = params
    Hkc = get_Hkc(Reθ, Hk, params, bltype)
    cτeq = Hs*(Hk - 1)*Hkc^2 / (2*Ga^2*Gb*(1 - Us)*H*Hk^2)
    return cτeq
end

delta(δs, θ, Hk) = min(3.15 + 1.72/(Hk - 1) + δs, 12*θ)  # TODO


# ----------------------------


# ---- averaging and upwinding ------
# factor used in upwinding
function upwind_factor(Hk1, Hk2, bltype)
    Cup = bltype == "wake" ? 5.0 : 1.0
    fhu = (Hk2 - 1) / (Hk1 - 1)
    eta = 1.0 - 0.5*exp(-log(abs(fhu))^2 * Cup / Hk2^2)
    return eta
end

# upwind two variables using eta factor
upwind(one, two, eta) = @. (1 - eta)*one + eta*two

# average two variables
average(one, two) = upwind(one, two, 0.5)


# -------- residuals ----------------

function lresid(state1, state2, x1, x2, params, bltype)

    # states
    θ1, δs1, v1, ue1 = state1
    θ2, δs2, v2, ue2 = state2
    θm, δsm, vm, uem = average(state1, state2)

    # positions
    xm = average(x1, x2)
    dx = x2 - x1

    # set third state
    if bltype == "laminar"
        n1 = v1; n2 = v2
        sqcτ1 = 1.0; sqcτ2 = 1.0  # unused
        cτ1 = 1.0; cτ2 = 1.0  # unused
    else
        n1 = 1.0; n2 = 1.0  # unused
        sqcτ1 = v1; sqcτ2 = v2
        cτ1 = sqcτ1^2; cτ2 = sqcτ2^2
    end

    # compressibility
    u1, M1, Reθ1 = compressibility(ue1, θ1, params)
    u2, M2, Reθ2 = compressibility(ue2, θ2, params)
    um, Mm, Reθm = compressibility(uem, θm, params)

    # shape factors
    H1 = δs1/θ1; H2 = δs2/θ2; Hm = δsm/θm

    Hk1, Hss1, Hs1 = shape(H1, M1, Reθ1, bltype)
    Hk2, Hss2, Hs2 = shape(H2, M2, Reθ2, bltype)
    Hkm, Hssm, Hsm = shape(Hm, Mm, Reθm, bltype)

    # upwind factor
    eta = upwind_factor(Hk1, Hk2, bltype)

    # skin friction
    cf1 = skinfriction(Reθ1, Hk1, M1, params.γ, bltype)
    cf2 = skinfriction(Reθ2, Hk2, M2, params.γ, bltype)
    cfm = skinfriction(Reθm, Hkm, Mm, params.γ, bltype)

    cfxt1 = cf1 * x1 / θ1
    cfxt2 = cf2 * x2 / θ2
    cfxtm = cfm * xm / θm
    cfxt = 0.25*cfxt1 + 0.25*cfxt2 + 0.5*cfxtm

    # momentum residual
    H = average(H1, H2)
    M = average(M1, M2)
    r_mom = log(θ2/θ1) + (2 + H - M^2)*log(u2/u1) - 0.5*log(x2/x1)*cfxt

    # dissipation
    Us1, cdi1 = dissipation(cτ1, Reθ1, H1, Hk1, Hs1, cf1, params.Gb, bltype)
    Us2, cdi2 = dissipation(cτ2, Reθ2, H2, Hk2, Hs2, cf2, params.Gb, bltype)

    # shape residual
    cfxtup = upwind(cfxt1, cfxt2, eta)
    cdxtHup = upwind((cdi1*x1)/(θ1*Hs1), (cdi2*x2)/(θ2*Hs2), eta)
    Hs = average(Hs1, Hs2)
    Hss = average(Hss1, Hss2)
    r_shape = log(Hs2/Hs1) + (2*Hss/Hs + 1 - H)*log(u2/u1) + log(x2/x1)*(cfxtup/2 - cdxtHup)


    # amplification residual (e^n method)
    if bltype == "laminar"

        dndx1 = amplification(θ1, n1, Reθ1, Hk1, params)
        dndx2 = amplification(θ2, n2, Reθ2, Hk2, params)
        dndx = average(dndx1, dndx2)

        r_amp = n2 - n1 - dndx*dx

        return [r_mom, r_shape, r_amp]

    else  # lag residual

        sqcτ = upwind(sqcτ1, sqcτ2, eta)
        cτeq1 = get_cτeq(Reθ1, H1, Hk1, Hs1, Us1)
        cτeq2 = get_cτeq(Reθ2, H2, Hk2, Hs2, Us2)
        sqcτeq = upwind(sqrt(cτeq1), sqrt(cτeq2), eta)

        δ = average(delta(δs1, θ1, Hk1), delta(δs2, θ2, Hk2))
        Us = average(Us1, Us2)

        cfup = upwind(cf1, cf2, eta)
        Hkup = upwind(Hk1, Hk2, eta)
        uq = get_uq(δsm, Reθm, Hkup, cfup, params, bltype)

        (; Klag, Gb, ηd) = params

        r_lag = 2*δ*log(sqcτ2/sqcτ1) - Klag/(Gb*(1 + Us))*(sqcτeq - ηd*sqcτ)*dx -
            2*δ*(uq*dx - log(u2/u1))

        return [r_mom, r_shape, r_lag]
    end

end