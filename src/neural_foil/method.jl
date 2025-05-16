#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    NeuralFoil <: Method

# Fields:
- `mean_inputs_scaled::Vector`
- `cov_inputs_scaled::Matrix`
- `inv_cov_inputs_scaled::Matrix`
- `weights::Vector{Vector}}`
- `biases::Vector{Vector}}`
- `Re::Float`
- `Ma::Float`
- `n_crit::Float'
- `xtr_upper::Float'
- `xtr_lower::Float'
"""
struct NeuralFoil{Tb,TM,TMa,Tn,TRe,TV,TW,Txl,Txu} <: Method
    mean_inputs_scaled::TV
    cov_inputs_scaled::TM
    inv_cov_inputs_scaled::TM
    weights::TW
    biases::Tb
    Re::TRe
    Ma::TMa
    n_crit::Tn
    xtr_upper::Txu
    xtr_lower::Txl
end

"""
    NeuralFoil(
        reynolds=1e6, mach=0.0; model_size="xlarge", n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0
    )

    Constructor for NeuralFoil type.

# Default Arguments:
- `reynolds::Float=1e6` : Reynolds number
- `mach::Float=0.0` : Mach number

# Keyword Arguments
- `model_size::String="xlarge"` : model size from NeuralFoil
- `n_crit::Float=9.0` : n_crit for Xfoil
- `xtr_upperFloat=1.0` : Xtr_Upper for Xfoil
- `xtr_lowerFloat=1.0` : Xtr_Lower for Xfoil

# Returns
- `method::NeuralFoil` : NeuralFoil method object
"""
function NeuralFoil(
    reynolds=1e6, mach=0.0; model_size="xlarge", n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0
)
    scaled_input_distribution = NPZ.npzread(
        joinpath(@__DIR__, "..", "..", "data", "scaled_input_distribution.npz")
    )
    nn_params = NPZ.npzread(
        joinpath(@__DIR__, "..", "..", "data", "nn-" * model_size * ".npz")
    )

    unique_keys = unique(sort([parse(Int, split(key, ".")[2]) for key in keys(nn_params)]))

    return NeuralFoil(
        scaled_input_distribution["mean_inputs_scaled"],
        scaled_input_distribution["cov_inputs_scaled"],
        scaled_input_distribution["inv_cov_inputs_scaled"],
        [nn_params["net.$(id).weight"] for id in unique_keys],
        [nn_params["net.$(id).bias"] for id in unique_keys],
        reynolds,
        mach,
        n_crit,
        xtr_upper,
        xtr_lower,
    )
end

"""
    NeuralOutputs

# Fields
- cl::Vector : lift coefficients
- cd::Vector : drag coefficients
- cm::Vector : moment coefficients
- confidence::Vector : confidence factor reported by NeuralFoil
"""
@kwdef struct NeuralOutputs{TM,TV}
    cl::TV
    cd::TV
    cm::TV
    confidence::TV
    top_xtr::TV
    bot_xtr::TV
    upper_bl_ue_over_vinf::TM
    lower_bl_ue_over_vinf::TM
    upper_theta::TM
    upper_H::TM
    lower_theta::TM
    lower_H::TM
end

#---------------------------------#
#            Neural Net           #
#---------------------------------#

"""
    squared_mahalanobis_distance(x::AbstractMatrix)

Computes the squared Mahalanobis distance of a set of points from the training data.

# Arguments:
- `x::Matrix`: Query point in the input latent space. Shape: (N_cases, N_inputs) For non-vectorized queries, N_cases=1.
- `mean_inputs_scaled::Vector` : from scaled_input_distribution
- `inv_cov_inputs_scaled::Matrix` : from scaled_input_distribution

# Returns:
- `sqmd::Vector` : The squared Mahalanobis distance. Shape: (N_cases,)
"""
function squared_mahalanobis_distance(x, mean_inputs_scaled, inv_cov_inputs_scaled)
    x_minus_mean = (x .- mean_inputs_scaled)'
    return sum((x_minus_mean * inv_cov_inputs_scaled) .* x_minus_mean; dims=2)
end

function sigmoid(x; ln_eps=log(10.0 / floatmax(Float64)))
    x_clipped = clamp.(x, ln_eps, -ln_eps)
    return 1.0 ./ (1.0 .+ exp.(-x_clipped))
end

function swish(x)
    return x ./ (1.0 .+ exp.(-x))
end

function net(x, cache)

    # Rename for Convenience
    (; weights, biases) = cache

    for (i, (W, b)) in enumerate(zip(weights, biases))
        x = W * x
        x .+= b
        if i != length(weights)
            x = swish(x)
        end
    end

    return x
end

#---------------------------------#
#         Mach Correction         #
#---------------------------------#

function smooth_beta(mach; blend_range=0.02)
    b(M) = 1.0 / sqrt(1.0 - min(M, 0.999)^2)
    return FLOWMath.quintic_blend(b(mach), b(0.99), mach, 0.975, blend_range)
end

function laitone_compressibility_correction(coeff, mach; gamma=1.4)
    beta = smooth_beta(mach)
    denom = beta + mach^2 * coeff / (2.0 * beta) * (1.0 + (gamma - 1) * mach^2 / 2.0)
    return coeff / denom
end

#---------------------------------#
#             Analysis            #
#---------------------------------#

function analyze_nf(coordinates, flow_angles; method=NeuralFoil())

    # Get CST parameters from coordinates
    cst_upper, cst_lower, cst_LE, cst_TE = at.determine_neuralfoil_cst(coordinates)

    # yu = at.neuralfoil_half_cst(cst_upper,xu,cst_TE/2,cst_LE)
    # yl = at.neuralfoil_half_cst(cst_lower,reverse(xl),-cst_TE/2,cst_LE)

    # Assemble inputs to run everything all at once
    x = stack([
        [
            cst_upper
            cst_lower
            cst_LE
            cst_TE * 50.0
            sind(2.0 * aoa)
            cosd(aoa)
            1.0 - cosd(aoa)^2
            (log(method.Re) - 12.5) / 3.5
            (method.n_crit .- 9.0) / 4.5
            method.xtr_upper
            method.xtr_lower
        ] for aoa in flow_angles
    ])

    x_flipped = stack([
        [
            -cst_lower
            -cst_upper
            -cst_LE
            cst_TE * 50.0
            -sind(2.0 * aoa)
            cosd(aoa)
            1.0 - cosd(aoa)^2
            (log(method.Re) - 12.5) / 3.5
            (method.n_crit .- 9.0) / 4.5
            method.xtr_lower
            method.xtr_upper
        ] for aoa in flow_angles
    ])

    # rename for convenience
    Wb = (; weights=method.weights, biases=method.biases)

    # - Call network (outputs num cases x num outputs) - #
    y = net(x, Wb)
    y_flipped = net(x_flipped, Wb)

    # - Compute confidence values - #
    y[1, :] .-=
        squared_mahalanobis_distance(
            x, method.mean_inputs_scaled, method.inv_cov_inputs_scaled
        ) ./ (2.0 * size(x, 1))
    y_flipped[1, :] .-=
        squared_mahalanobis_distance(
            x_flipped, method.mean_inputs_scaled, method.inv_cov_inputs_scaled
        ) ./ (2.0 * size(x_flipped, 1))

    # - Unflip flipped output - #
    y_unflipped = copy(y_flipped)
    y_unflipped[2, :] .*= -1.0  # CL
    y_unflipped[4, :] .*= -1.0  # CM
    y_unflipped[5, :] .= y_flipped[6, :]   # Top_Xtr
    y_unflipped[6, :] .= y_flipped[5, :]   # Bot_Xtr

    # Switch upper and lower Ret, H
    y_unflipped[7:(7 + 32 * 2 - 1), :] .= y_flipped[(7 + 32 * 3):(7 + 32 * 5 - 1), :]
    y_unflipped[(7 + 32 * 3):(7 + 32 * 5 - 1), :] .= y_flipped[7:(7 + 32 * 2 - 1), :]

    # Switch upper_bl_ue/vinf with lower_bl_ue/vinf
    y_unflipped[(7 + 32 * 2):(7 + 32 * 3 - 1), :] .=
        -y_flipped[(7 + 32 * 5):(7 + 32 * 6 - 1), :]
    y_unflipped[(7 + 32 * 5):(7 + 32 * 6 - 1), :] .=
        -y_flipped[(7 + 32 * 2):(7 + 32 * 3 - 1), :]

    # - Average outputs - #
    y_fused = (y .+ y_unflipped) ./ 2.0
    y_fused[1, :] .= sigmoid.(y_fused[1, :])
    y_fused[5, :] .= clamp.(y_fused[5, :], 0, 1)
    y_fused[6, :] .= clamp.(y_fused[6, :], 0, 1)

    # Set up outputs for return
    N = 32 #hard coded in neuralfoil
    confidence = y_fused[1, :]
    cl = y_fused[2, :] ./ 2.0
    cd = exp.((y_fused[3, :] .- 2.0) .* 2)
    cm = y_fused[4, :] ./ 20.0
    top_xtr = y_fused[5, :]
    bot_xtr = y_fused[6, :]
    upper_bl_ue_over_vinf = y_fused[(7 + N * 2):(7 + N * 3 - 1), :]
    lower_bl_ue_over_vinf = y_fused[(7 + N * 5):(7 + N * 6 - 1), :]
    upper_theta =
        ((10.0 .^ y_fused[7:(7 + N - 1), :]) .- 0.1) ./
        (abs.(upper_bl_ue_over_vinf) .* method.Re)
    upper_H = 2.6 .* exp.(y_fused[(7 + N):(7 + N * 2 - 1), :])
    lower_theta =
        ((10.0 .^ y_fused[(7 + N * 3):(7 + N * 4 - 1), :]) .- 0.1) ./
        (abs.(lower_bl_ue_over_vinf) .* method.Re)
    lower_H = 2.6 .* exp.(y_fused[(7 + N * 4):(7 + N * 5 - 1), :])

    # - Apply Mach Corrections - #
    if !iszero(method.Ma)
        # Laitone Compressibility correction
        cl = laitone_compressibility_correction.(cl, method.Ma)
        cm = laitone_compressibility_correction.(cm, method.Ma)
    end

    # - Return Outputs - #
    return NeuralOutputs(
        cl,
        cd,
        cm,
        confidence,
        top_xtr,
        bot_xtr,
        upper_bl_ue_over_vinf,
        lower_bl_ue_over_vinf,
        upper_theta,
        upper_H,
        lower_theta,
        lower_H,
    )
end
