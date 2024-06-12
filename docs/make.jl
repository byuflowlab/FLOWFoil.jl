using Documenter
using FLOWFoil

makedocs(;
    modules=[FLOWFoil, FLOWFoil.AirfoilTools],
    authors="Judd Mehr",
    repo="https://github.com/byuflowlab/FLOWFoil.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        repolink="https://github.com/byuflowlab/FLOWFoil.jl/blob/{commit}{path}#L{line}",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "FLOWFoil" => [
            "Getting Started" => "FLOWFoil/tutorial.md",
            "Examples" => "FLOWFoil/examples.md",
            "API Index" => "FLOWFoil/reference.md",
            "Theory" => "FLOWFoil/theory.md",
        ],
        "AirfoilTools" => [
            "Intro" => "AirfoilTools/intro.md",
            "Airfoil Generation" => "AirfoilTools/parameterizations.md",
            "Airfoil Manipulation" => "AirfoilTools/geometry_manipulations.md",
            "API Reference" => "AirfoilTools/api.md",
        ],
    ],
    sitename="FLOWFoil.jl",
    authors="Judd Mehr <juddmehr@byu.edu>",
    checkdocs=:exports,
)

deploydocs(; repo="github.com/byuflowlab/FLOWFoil.jl", devbranch="dev")
