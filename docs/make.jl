using FLOWFoil
using Documenter

DocMeta.setdocmeta!(FLOWFoil, :DocTestSetup, :(using FLOWFoil); recursive=true)

makedocs(;
    modules=[FLOWFoil],
    authors="Judd Mehr",
    repo="https://github.com/byuflowlab/FLOWFoil.jl/blob/{commit}{path}#{line}",
    sitename="FLOWFoil.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://flow.byu.edu/FLOWFoil.jl",
        assets=String[],
    ),
    pages=[
        "Intro" => "index.md",
        "Quick Start" => "tutorial.md",
        "Guided Examples" => "examples.md",
        "API Reference" => "reference.md",
        "Theory" => "theory.md",
    ],
)

deploydocs(; repo="github.com/byuflowlab/FLOWFoil.jl", devbranch="main")
