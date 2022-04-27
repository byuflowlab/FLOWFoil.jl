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
        canonical="https://byuflowlab.github.io/FLOWFoil.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/FLOWFoil.jl",
    devbranch="main",
)
