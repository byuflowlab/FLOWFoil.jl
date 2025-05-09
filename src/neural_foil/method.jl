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

function analyze_nf(
    coordinates, flow_angles; reynolds=[0.0], machs=[0.0], method=NeuralFoil()
)
    aero = get_aero_from_coordinates(
        np_array(reverse(coordinates; dims=1));
        alpha=flow_angles,
        Re=reynolds,
        model_size=method.model_size,
    )
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
end
