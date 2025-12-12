# Calling DuctAPE from Python

In this example we repeat the [quick start](tutorial.md), but with Python.  We won't repeat all the usage details described in the quick start.  We are using the [JuliaCall](https://juliapy.github.io/PythonCall.jl/stable/) packages to enable this.

## Setup

```python
# optional: point python to a specific julia installation
# import os
# os.environ["JULIA_BINDIR"] = path_to_julia

from juliacall import Main as jl

# If needed, install DuctAPE into the python julia environment
jl.seval("using Pkg")
jl.Pkg.add("FLOWFoil")

# Load DuctAPE
jl.seval("using FLOWFoil")

import numpy as np
```

## Example

```python
# - Same Quick Start as Before, but with some minor changes for Python syntax - #

# AirfoilTools is exported from FLOWFoil
x, y = jl.AirfoilTools.naca4()

# choose angles of attack
angles_of_attack = np.arange(-5.0, 15.0 + 1e-6, 1.0)  # include 15.0

# Run solve
outputs = jl.analyze(x, y, angles_of_attack)

# Print some outputs.
println("cl: ", outputs.cl)
```

## Automatic Derivatives

We now continue the example demonstrating how to get derivatives for use in Python (but where the derivative computation occurs in Julia via algorithmic differentiation).  For this functionality we need to load the [PythonCall](https://juliapy.github.io/PythonCall.jl/stable/) package (which enables us to call back into Python from Julia), and we need the [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl) packages which provides a convenience function for the derivative computation.  Alternatively, for more advanced users you can just use the Julia differentiation packages directly (ForwardDiff, ReverseDiff, etc.).  Note that for Julia AD to work, all the function calls will need to be Julia function calls.  Even though all the setup is happening here in Python, we are only setting up inputs.  All functions are calls to Julia (jl.somefunction())

The below example demonstrates a forward mode Jacobian and a forward mode Jacobian-vector product.  Other options (reverse Jacobian, reverse vector-Jacobian product) are discussed in [ImplicitAD](https://github.com/byuflowlab/ImplicitAD.jl) docs.

```python
# load PythonCall BEFORE ImplicitAD
jl.seval("using PythonCall")

# Install ImplicitAD as needed
jl.seval("using Pkg")
jl.Pkg.add("ImplicitAD")

# Load ImplicitAD.derivativesetup
jl.seval("using ImplicitAD: derivativesetup")

nc = len(x)

# ImplicitAD expects a funciton we want to differentiate in the form f = func(x, p)
# where f is output vector, x is input vector, and p are parameters we do not differentiate w.r.t.
def ffwrap(x, p):

    outputs = jl.analyze(x, y, angles_of_attack)

    return [outputs.cl]


p = ()
jacobian = jl.derivativesetup(
    ffwrap, x, p, "fjacobian"
)  # a forward-mode Jacobian is one option

# preallocate Jacobian then evaluate
J = np.zeros((2, len(x)))
jacobian(J, x)
print(J)
# can now change x, and evaluate jacobian(J, x) repeatedly at other points

# demonstrate a Jacobian-vector product
jvp = jl.derivativesetup(ffwrap, x, p, "jvp")
xdot = np.ones(len(x))
fdot = np.zeros(2)
jvp(fdot, x, xdot)
print(fdot)
# can continue to call jvp for different x, xdot pairs
```
