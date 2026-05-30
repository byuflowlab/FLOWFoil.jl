"""
    Xfoil <: Method

Configuration for an Xfoil-like panel method. When `viscous=false` (default), only the
inviscid linear-vortex panel solve is run. When `viscous=true`, the coupled viscous
boundary-layer + inviscid system is solved (single-body only; multi-alpha is supported
by looping internally).

Viscous-only fields are ignored when `viscous=false`.

# Fields
- `viscous::Bool` : run inviscid (`false`) or coupled viscous (`true`) analysis.

# Viscous-only fields
- `Re` : Reynolds number based on chord and freestream.
- `Mach` : freestream Mach number.
- `rho` : non-dimensional freestream density.
- `Vinf` : non-dimensional freestream speed.
- `chord` : reference chord length.
- `ncrit` : critical amplification factor for free transition (e^n method).
- `wake_length` : wake length in chords.
- `numpanels_wake` : number of wake panels (set automatically from airfoil node count if `nothing`).
- `etol` : Newton convergence tolerance.
- `maxiters` : maximum Newton iterations (both local station solves and the coupled global solve).
- `xft_l` : forced lower-surface transition x/c (set ≥ 1 to disable forcing).
- `xft_u` : forced upper-surface transition x/c (set ≥ 1 to disable forcing).
- `verbose` : print iteration info during viscous solve.
"""
@kwdef struct Xfoil{TB,TF,TI,TIW} <: Method
    viscous::TB = false
    Re::TF = 1e6
    Mach::TF = 0.0
    rho::TF = 1.0
    Vinf::TF = 1.0
    chord::TF = 1.0
    ncrit::TF = 9.0
    wake_length::TF = 1.0
    numpanels_wake::TIW = nothing
    etol::TF = 1e-10
    maxiters::TI = 50
    xft_l::TF = 1.0
    xft_u::TF = 1.0
    verbose::TB = false
end

const Mfoil = Xfoil
