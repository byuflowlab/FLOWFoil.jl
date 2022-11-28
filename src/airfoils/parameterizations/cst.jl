#=
Class-Shape Transformation (CST) Airfoil Parameterization

Authors: Judd Mehr,

=#


function bernstein(r, n, x)
    return binomial(n, r) .* x.^r .* (1 .- x).^(n.-r)
end

function cst(coeffs::AbstractArray{<:Number,1}, N::Integer=80,
    x::AbstractArray{<:Number,1}=cosinespacing(N), dzu::Number=0.0, dzl::Number=0.0,
    N1::Number=0.5, N2::Number=1.0)

    if length(coeffs)%2 != 0
        error("CST: Must have even number of coefficients (half upper, half lower)")
    else
        n = Int(length(coeffs)/2)
    end

    coeffU = coeffs[1:n]
    coeffL = coeffs[n+1:end]

    zu = halfcst(coeffU, x, dzu, N1, N2)
    zl = halfcst(coeffL, x, dzl, N1, N2)

    return x, zu, zl

end

function halfcst(coeffs::AbstractArray{<:Number,1}, x::AbstractArray{<:Number,1}=cosinespacing(N),
    dz::Number=0.0, N1::Number=0.5, N2::Number=1.0)

    n = length(coeffs)

    C = x.^N1.*(1 .- x).^N2

    nb = length(coeffs)-1

    nx = length(x)

    S = zeros(eltype(coeffs), nx)

    for i=1:nb+1
        S += coeffs[i]*bernstein(i-1, nb, x)
    end

    z = C .* S + x*dz

    return z

end

function getcst(x::AbstractArray{<:Number, 1}, zu::AbstractArray{<:Number, 1},
    zl::AbstractArray{<:Number, 1}, n::Integer, dzu::Number=0.0, dzl::Number=0.0,
    N1::Number=0.5, N2::Number=1.0)

    # models to fit
    cstU(x, coeffU) = halfcst(coeffU, x, dzu, N1, N2)
    cstL(x, coeffL) = halfcst(coeffL, x, dzl, N1, N2)

    # initial guesses
    coeffU = ones(n)
    coeffL = -ones(n)

    # solve for coefficients
    ufit = LsqFit.curve_fit(cstU, x, zu, coeffU)
    lfit = LsqFit.curve_fit(cstL, x, zl, coeffL)

    return vcat(ufit.param, lfit.param)
end
