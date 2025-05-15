"""
    NeuralFoil <: Method

# Fields:
- model_size::String = "xlarge" : the "model size" for NeuralFoil. Choose from: "xxsmall", "xsmall", "small", "medium", "large", "xlarge", "xxlarge", and "xxxlarge".
"""
@kwdef struct NeuralFoil{TS} <: Method
    model_size::TS = "xlarge"
end

"""
    NeuralOutputs

# Fields
- cl::Vector : lift coefficients
- cd::Vector : drag coefficients
- cm::Vector : moment coefficients
- confidence::Vector : confidence factor reported by NeuralFoil
"""
@kwdef struct NeuralOutputs{TV}
    cl::TV
    cd::TV
    cm::TV
    confidence::TV
end

function smooth_beta(mach; blend_range=0.02)
    b(M) = 1.0 / sqrt(1.0 - min(M, 0.999)^2)

    return FLOWMath.quintic_blend(b(mach), b(0.99), mach, 0.975, blend_range)
end

function laitone_compressibility_correction(coeff, mach; gamma=1.4)
    beta = smooth_beta(mach)

    denom = beta + mach^2 * coeff / (2.0 * beta) * (1.0 + (gamma - 1) * mach^2 / 2.0)

    return coeff / denom
end

function nf_val(variable_inputs, constant_parameters)
    pyx = pyconvert(Py, variable_inputs)
    pyp = pyconvert(Py, constant_parameters)

    return pyconvert(Vector{Float64}, nf_val_py_wrap(pyx, pyp; call=true))
end

function nf_jvp(variable_inputs, constant_parameters, push_vector)
    pyx = pyconvert(Py, variable_inputs)
    pyv = pyconvert(Py, push_vector)
    pyp = pyconvert(Py, constant_parameters)

    return pyconvert(Vector{Float64}, nf_jvp_py_wrap(pyx, pyv, pyp))
end

function nf_vjp(variable_inputs, constant_parameters, pull_covector)
    pyx = pyconvert(Py, variable_inputs)
    pyc = pyconvert(Py, pull_covector)
    pyp = pyconvert(Py, constant_parameters)

    return pyconvert(Vector{Float64}, nf_vjp_py_wrap(pyx, pyc, pyp))
end

function analyze_nf(coordinates, flow_angles; reynolds=0.0, mach=0.0, method=NeuralFoil())
    # - Pack Inputs - #

    # vectorize inputs
    variable_inputs = [reduce(vcat, reverse(coordinates; dims=1)); flow_angles; reynolds]

    # pass in sizing for unpacking in JAX
    lfa = length(flow_angles)
    constant_parameters = (;
        coord_length=length(coordinates),
        coord_shape=size(coordinates),
        angle_length=lfa,
        model_size=method.model_size,
    )

    # - Run through ImplicitAD - #
    output_vector_with_duals = ImplicitAD.provide_rule(
        nf_val, variable_inputs, constant_parameters; mode="vp", jvp=nf_jvp, vjp=nf_vjp
    )

    # - Unpack Outputs - #
    cl = output_vector_with_duals[1:lfa]
    cd = output_vector_with_duals[(lfa + 1):(lfa * 2)]
    cm = output_vector_with_duals[(lfa * 2 + 1):(lfa * 3)]
    confidence = output_vector_with_duals[(lfa * 3 + 1):(lfa * 4)]

    # - Apply Mach Corrections - #
    if iszero(mach)
        return NeuralOutputs(; cl=cl, cd=cd, cm=cm, confidence=confidence)
    else

        # Apply Mach Corrections if Mach != 0.0
        # Laitone Compressibility correction
        cl = laitone_compressibility_correction.(cl, mach)
        cm = laitone_compressibility_correction.(cm, mach)
        return NeuralOutputs(; cl=cl, cd=cd, cm=cm, confidence=confidence)
    end
end
