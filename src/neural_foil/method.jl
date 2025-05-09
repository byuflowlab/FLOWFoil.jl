@kwdef struct NeuralFoil{TS} <: Method
    model_size::TS = "xlarge"
end

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
        cl=aero["CL"], cd=aero["CD"], cm=aero["CM"], confidence=aero["analysis_confidence"]
    )
end
