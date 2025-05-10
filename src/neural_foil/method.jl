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

function analyze_nf(coordinates, flow_angles; reynolds=0.0, mach=0.0, method=NeuralFoil())

    #TODO: write a wrapper function for jax.jit usage that takes in an array of inputs and returns an array of outputs.
    #TODO: figure out how to have those inputs contain derivative information
    #TODO: figure out how to have those outputs contain derivative information as well
    #    for outputs, see: https://docs.jax.dev/en/latest/notebooks/autodiff_cookbook.html#evaluate-a-function-and-its-gradient-using-value-and-grad
    #    or rather: https://docs.jax.dev/en/latest/notebooks/autodiff_cookbook.html#jacobians-and-hessians-using-jacfwd-and-jacrev
    #    or probably (if using implicitAD) https://docs.jax.dev/en/latest/notebooks/autodiff_cookbook.html#jvps-in-jax-code
    #    general page: https://docs.jax.dev/en/latest/notebooks/autodiff_cookbook.html#jvps-in-jax-code
    #TODO: does the whole function need to be in python?  if so, probably need to add a python module with all the functions needed and then figure out how to load that manually with the CondaPkg stuff in order to PythonCall macro the function call into julia

    aero = get_aero_from_coordinates(
        jnp_array(reverse(coordinates; dims=1));
        alpha=jnp_array(flow_angles),
        Re=reynolds,
        model_size=method.model_size,
    )

    if iszero(mach)
        return NeuralOutputs(;
            cl=pyconvert(Array{Float64}, aero["CL"]),
            cd=pyconvert(Array{Float64}, aero["CD"]),
            cm=pyconvert(Array{Float64}, aero["CM"]),
            confidence=pyconvert(Array{Float64}, aero["analysis_confidence"]),
            # cl=aero["CL"],
            # cd=aero["CD"],
            # cm=aero["CM"],
            # confidence=aero["analysis_confidence"],
        )
    else

        # Apply Mach Corrections if Mach != 0.0
        # Laitone Compressibility correction
        cl =
            laitone_compressibility_correction.(pyconvert(Array{Float64}, aero["CL"]), mach)
        cm =
            laitone_compressibility_correction.(pyconvert(Array{Float64}, aero["CM"]), mach)

        return NeuralOutputs(;
            cl=cl,
            cd=pyconvert(Array{Float64}, aero["CD"]),
            cm=cm,
            confidence=pyconvert(Array{Float64}, aero["analysis_confidence"]),
            # cl=aero["CL"],
            # cd=aero["CD"],
            # cm=aero["CM"],
            # confidence=aero["analysis_confidence"],
        )
    end
end
