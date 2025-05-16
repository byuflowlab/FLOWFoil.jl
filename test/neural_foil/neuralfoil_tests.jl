using PythonCall
using CondaPkg
CondaPkg.add_pip("neuralfoil")
CondaPkg.add_pip("aerosandbox")

__precompile__(false)

macro wrap_pyfuns(modsym, fname, cname)
    quote
        const pymod = pyimport($modsym)
        # Define functions with the same names as the input symbols
        $(:(function $(esc(cname))(args...; kwargs...)
            pyf = @pyconst pymod.$(fname)
            return pyf(args...; kwargs...)
        end))
    end
end

@wrap_pyfuns "neuralfoil" get_aero_from_coordinates get_aero_from_coordinates
@wrap_pyfuns "aerosandbox.geometry.airfoil.airfoil_families" get_kulfan_parameters get_kulfan_parameters
@wrap_pyfuns "numpy" array np_array

function wrapped_nf(
    coordinates, flow_angles; reynolds=[0.0], machs=[0.0], model_size="xlarge"
)
    aero = get_aero_from_coordinates(
        np_array(reverse(coordinates; dims=1));
        alpha=flow_angles,
        Re=reynolds,
        model_size=model_size,
    )
    return (;
        cl=pyconvert(Vector{Float64}, aero["CL"]),
        cd=pyconvert(Vector{Float64}, aero["CD"]),
        cm=pyconvert(Vector{Float64}, aero["CM"]),
        confidence=pyconvert(Vector{Float64}, aero["analysis_confidence"]),
    )
end

@testset "Compare Kulfan Solvers" begin
    x, y = at.naca4(2.0, 4.0, 12.0)
    coordinates = [x y]

    cst_outs = get_kulfan_parameters(np_array(reverse(coordinates; dims=1)))
    cst_lower_nf = pyconvert(Vector{Float64}, cst_outs["lower_weights"])
    cst_upper_nf = pyconvert(Vector{Float64}, cst_outs["upper_weights"])
    cst_TE_nf = pyconvert(Float64, cst_outs["TE_thickness"])
    cst_LE_nf = pyconvert(Float64, cst_outs["leading_edge_weight"])

    cst_upper, cst_lower, cst_LE, cst_TE = at.determine_neuralfoil_cst(coordinates)

    @test isapprox(cst_upper, cst_upper_nf, atol=1e-8)
    @test isapprox(cst_lower, cst_lower_nf, atol=1e-8)
    @test isapprox(cst_TE, cst_TE_nf, atol=1e-8)
    @test isapprox(cst_LE, cst_LE_nf, atol=1e-8)
end

@testset "Compare to NeuralFoil" begin
    x, y = at.naca4(2.0, 4.0, 12.0)
    coordinates = [x y]
    flow_angles = range(-5, 15; step=1)
    reynolds = 1e6
    mach = 0.0
    model_size = "xlarge"

    outputs_nfpy = wrapped_nf(
        coordinates, flow_angles; reynolds=reynolds, machs=mach, model_size=model_size
    )

    outputs_ff = analyze(
        coordinates, flow_angles; method=NeuralFoil(reynolds, mach; model_size=model_size)
    )

    @test isapprox(outputs_nfpy.cl, outputs_ff.cl, atol=1e-8)
    @test isapprox(outputs_nfpy.cd, outputs_ff.cd, atol=1e-8)
    @test isapprox(outputs_nfpy.cm, outputs_ff.cm, atol=1e-8)
end
