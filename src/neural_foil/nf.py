import jax
import jax.numpy as jnp
import neuralfoil as nf


def nf_val_py(variable_inputs, constant_parameters, call=False):
    """
    Run NeuralFoil.get_aero_from_coordinates.
    """

    def f(variable_inputs):

        # - Unpack Variables - #
        variable_inputs = jnp.array(variable_inputs)
        coordinates = jnp.column_stack(
            (
                variable_inputs[0 : constant_parameters.coord_shape[0]],
                variable_inputs[
                    constant_parameters.coord_shape[0] : 2
                    * constant_parameters.coord_shape[0]
                ],
            )
        )

        flow_angles = variable_inputs[
            constant_parameters.coord_length : constant_parameters.coord_length
            + constant_parameters.angle_length
        ]

        reynolds = variable_inputs[-1]

        # - Run NeuralFoil - #
        aero = nf.get_aero_from_coordinates(
            coordinates,
            alpha=flow_angles,
            Re=reynolds,
            model_size=constant_parameters.model_size,
        )

        # - Return concatenated results - #
        return jnp.concatenate(
            [
                jnp.atleast_1d(aero["CL"]),
                jnp.atleast_1d(aero["CD"]),
                jnp.atleast_1d(aero["CM"]),
                jnp.atleast_1d(aero["analysis_confidence"]),
            ]
        )

    if call:
        print("call true")
        return f(variable_inputs)
    else:
        print("call false")
        return f


def nf_jvp_py(variable_inputs, push_vector, constant_parameters):
    print("jvp")
    x = jnp.array(variable_inputs)
    v = jnp.array(push_vector)
    f = nf_val_py(x, constant_parameters)
    y, jvp = jax.jvp(f, (x,), (v,))
    return jvp


def nf_vjp_py(variable_inputs, pull_covector, constant_parameters):
    print("jvp")
    x = jnp.array(variable_inputs)
    c = jnp.array(pull_covector)
    f = nf_val_py(x, constant_parameters)
    y, vjp_fun = jax.vjp(f, x)
    vjp = vjp_fun(c)[0]
    return vjp
