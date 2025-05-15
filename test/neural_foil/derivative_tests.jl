@testset "NeuralFoil Wrapper Derivatives" begin
    function vpfun(np)

        # generate geometry
        x, y = AirfoilTools.naca4(np[1], np[2], np[3])
        coordinates = [x y]

        flow_angles = range(-5.0, 5.0; step=5)
        reynolds = 1e6
        model_size = "xxxlarge"
        method = FLOWFoil.NeuralFoil(; model_size=model_size)

        # vectorize inputs
        variable_inputs = [
            reduce(vcat, reverse(coordinates; dims=1))
            flow_angles
            reynolds
        ]

        # pass in sizing for unpacking in JAX
        constant_parameters = (;
            coord_length=length(coordinates),
            coord_shape=size(coordinates),
            angle_length=length(flow_angles),
            model_size=method.model_size,
        )

        return ImplicitAD.provide_rule(
            FLOWFoil.nf_val,
            variable_inputs,
            constant_parameters;
            mode="vp",
            jvp=FLOWFoil.nf_jvp,
            vjp=FLOWFoil.nf_vjp,
        )
    end

    function fdfun(np)

        # generate geometry
        x, y = AirfoilTools.naca4(np[1], np[2], np[3])
        coordinates = [x y]

        flow_angles = range(-5.0, 5.0; step=5)
        reynolds = 1e6
        model_size = "xxxlarge"
        method = FLOWFoil.NeuralFoil(; model_size=model_size)

        # vectorize inputs
        variable_inputs = [
            reduce(vcat, reverse(coordinates; dims=1))
            flow_angles
            reynolds
        ]

        # pass in sizing for unpacking in JAX
        constant_parameters = (;
            coord_length=length(coordinates),
            coord_shape=size(coordinates),
            angle_length=length(flow_angles),
            model_size=method.model_size,
        )

        return ImplicitAD.provide_rule(
            FLOWFoil.nf_val, variable_inputs, constant_parameters; mode="cfd"
        )
    end

    vpjac = ForwardDiff.jacobian(vpfun, [2.0, 4.0, 12.0])
    fdjac = ForwardDiff.jacobian(fdfun, [2.0, 4.0, 12.0])

    @test isapprox(vpjac, fdjac, atol=1e-5)
end
