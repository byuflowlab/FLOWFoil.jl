#=
Closure equations

Authors: Judd Mehr,

Date Started: 23 May 2022

Change Log:
=#

# SHAPE PARAMETER CLOSURES
"""
"""
function Hk(H, Me; limit=false, wake=false)
    hk = (H - 0.29Me^2) / (1.0 + 0.113 * Me^2)
    if limit
        if wake
            return max(hk, 1.00005)
        else
            return max(hk, 1.05)
        end
    else
        return hk
    end
end

"""
"""
function Hstarstar(Me, Hk)
    return 0.064Me^2 / (Hk - 0.8)
end

"""
"""
function Hstarlam(Hk)
    Htildek = Hk - 4.35
    if Hk < 4.35
        return (0.0111 * Htildek^2 - 0.0278 * Htildek^3) / (Hk - 1.0) + 1.528 -
               0.002 * (Hk * Htildek)^2
    else
        return 0.015 * Htildek^2 / Hk + 1.528
    end
end

"""
"""
function Hstarturbinc(Retheta, Hk)
    H0 = min(3.0 + 400.0 / Retheta, 4.0)
    Hr = (H0 - Hk) / (H0 - 1.0)
    Retildetheta = max(Retheta, 200.0)
    AHstar = Hk - H0 + 4.0 / log(Retildetheta)
    if Hk < H0
        return 1.5 +
               4.0 / Retildetheta +
               1.5 * (0.5 - 4.0 / Retildetheta) * Hr^2.0 / (Hk + 0.5)
    else #TODO: unclear if AHstar^2 goes inside log.
        return 1.5 +
               4.0 / Retildetheta +
               (Hk - H0)^2 * (0.007 * log(Retildetheta) / AHstar^2 + 0.015 / Hk)
    end
end

"""
"""
function delta(Hk, deltastar, theta)
    return min(3.15 + 1.72 / (Hk - 1.0) + deltastar, 12.0 * theta)
end

"""
"""
function Hstarturb(Hstarturbinc, Me)
    return (Hstarturbinc + 0.028 * Me^2) / (1.0 + 0.014 * Me^2)
end

# SKIN FRICTION COEFFICIENT
"""
"""
function cflam(Retheta, Hk)
    if Hk < 5.5
        return (0.0727 * (5.5 - Hk)^3 / (Hk + 1.0) - 0.07) / Retheta
    else
        return (0.015 * (1.0 - 1.0 / (Hk - 4.5))^2 - 0.07) / Retheta
    end
end

"""
"""
function Acf(Hk)
    acf = -1.33 * Hk
    if acf < -17.0
        return -20.0 + 3.0 * exp((acf + 17.0) / 3.0)
    else
        return acf
    end
end

"""
"""
function Fc(Me, gamma_air)
    return sqrt(1.0 + 0.5 * (gamma_air - 1.0) * Me^2)
end

"""
"""
function Bcf(Retheta, Fc)
    return max(log(10, Retheta / Fc), 1.303)
end

"""
"""
function cfturb(Acf, Bcf, Hk)
    return 0.3 * exp(Acf) * Bcf^(-1.74 - 0.31 * Hk) +
           0.00011 * (tanh(4.0 - Hk / 0.875) - 1.0)
end

# DISSIPATION COEFFICIENT
"""
"""
function cdilam(Retheta, Hk)
    if Hk < 4.0
        return (0.00205 * (4.0 - Hk)^5.5 + 0.207) / Retheta
    else
        return (-0.0016 * (Hk - 4.0)^2 / (1.0 + 0.02 * (Hk - 4.0)^2) + 0.207) / Retheta
    end
end

"""
"""
function Us(Hstar, Hk, H, Gbeta; wake=false)
    us = 0.5 * Hstar * (1.0 - (Hk - 1.0) / (H * Gbeta))

    if wake
        return min(us, 0.99995)
    else
        return min(us, 0.98)
    end
end
"""
"""
function cdiwall(cf, Us, Hstar, Hk, Retheta)
    return 0.5 *
           cf *
           Us *
           (2.0 / Hstar) *
           0.5 *
           (1.0 + tanh((Hk - 1.0) * log(Retheta) / 2.1))
end

"""
"""
function cdiouter(ctau, Us, Hstar)
    return ctau * (0.995 - Us) * 2.0 / Hstar
end

"""
"""
function cdistress(Us, Hstar, Retheta)
    return 0.3 * (0.995 - Us)^2 / (Hstar * Retheta)
end

"""
"""
function cditurb(cdiwall, cdiouter, cdistress, cdilam)
    return min(cdiwall + cdiouter + cdistress, cdilam)
end

"""
"""
function cdilamwake(Hk, Hstar, Retheta)
    return 2.2 * (1.0 - 1.0 / Hk)^2 * (1.0 / Hk) / (Hstar * Retheta)
end

"""
"""
function cdiwake(cdiouter, cdistress, cdilamwake)
    return min(cdiouter + cdistress, cdilamwake)
end

# AMPLIFICATION FACTOR INCREASE RATE
"""
"""
function epsntilde(ntilde, ntildecrit)
    return 0.001 * (1.0 + tanh(5.0 * (ntilde - ntildecrit)))
end

"""
"""
function Hhat(Hk)
    return 1.0 / (Hk - 1.0)
end

"""
"""
function L0(Hhat)
    return 2.492 * Hhat^0.43 + 0.7 * (1.0 + tanh(14.0 * Hhat - 9.24))
end

"""
"""
function sntilde(Retheta, L0)
    return (log(10, Retheta) - (L0 - 0.1)) / 0.2
end

"""
"""
function rntilde(sntilde)
    if sntilde < 0.0
        return 0.0
    elseif sntilde > 1.0
        return 1.0
    else
        return 3.0 * sntilde^2 - 2.0 * sntilde^3
    end
end

"""
"""
function gntilde(Hhat)
    return 0.028 / Hhat - 0.0345 * exp((-3.87 * Hhat - 2.52)^2)
end

"""
"""
function fntilde(Hhat)
    return -0.05 + 2.7 * Hhat - 5.5 * Hhat^2 + 3.0 * Hhat^3 + 0.1 * exp(-20.0 * Hhat)
end

"""
"""
function dntildedxi(rntilde, fntilde, gntilde, epsntilde, theta)
    return (rntilde * fntilde * gntilde + epsntilde) / theta
end

# EQUILIBRIUM VALUES
"""
"""
function Hkc(Hk, Gc, Retheta; wake=false)
    Gtildec = wake ? 0.0 : Gc
    return Hk - 1.0 - Gtildec / Retheta
end

"""
"""
function uq(cf, Hkc, GA, etad, Hk, Gbeta, deltastar)
    return (0.5 * cf - (Hkc / (GA * etad * Hk))^2) / (Gbeta^deltastar)
end

#TODO: check mfoil code for typos...
"""
"""
function ctaueq(H, Hstar, Hk, Hkc, GA, Gbeta, Us)
    return Hstar * (Hk - 1.0) * Hkc^2 / (2.0 * GA^2 * Gbeta * (1.0 - Us) * H * Hk^2)
end
