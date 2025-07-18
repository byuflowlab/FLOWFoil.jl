using Documenter
using FLOWFoil

makedocs(;
    modules=[FLOWFoil, FLOWFoil.AirfoilTools],
    repo="https://github.com/byuflowlab/FLOWFoil.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        repolink="https://github.com/byuflowlab/FLOWFoil.jl/blob/{commit}{path}#L{line}",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "FLOWFoil" => [
            "Quick Start" => "FLOWFoil/tutorial.md",
            "Tutorials" => "FLOWFoil/additional_tutorials.md",
            "Additional Examples" => "FLOWFoil/examples.md",
            "API Index" => "FLOWFoil/api.md",
        ],
        "AirfoilTools" => [
            "Intro" => "AirfoilTools/intro.md",
            "Airfoil Generation" => "AirfoilTools/parameterizations.md",
            "Airfoil Manipulation" => "AirfoilTools/geometry_manipulations.md",
            "API Index" => "AirfoilTools/api.md",
        ],
    ],
    sitename="FLOWFoil.jl",
    authors="Judd Mehr <juddmehr@byu.edu>, Tyler Critchfield, Ayden Bennett <jbenne26@byu.edu>, Nathan Lehnhof",
    # checkdocs=:exports,
    checkdocs=:none,
)

deploydocs(; repo="github.com/byuflowlab/FLOWFoil.jl", devbranch="main")
