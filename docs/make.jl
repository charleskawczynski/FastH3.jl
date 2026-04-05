using Documenter
using FastH3

makedocs(;
    sitename="FastH3.jl",
    modules=[FastH3],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://charleskawczynski.github.io/FastH3.jl",
    ),
    clean=true,
    pages=[
        "Home" => "index.md",
        "API Reference" => [
            "Indexing" => "api/indexing.md",
            "Inspection" => "api/inspection.md",
            "Traversal" => "api/traversal.md",
            "Hierarchy" => "api/hierarchy.md",
            "Directed Edges" => "api/directed_edges.md",
            "Vertexes" => "api/vertexes.md",
            "Measurement" => "api/measurement.md",
        ],
        "Internals" => "api/internals.md",
    ],
    doctest=true,
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/charleskawczynski/FastH3.jl.git",
    target="build",
    devbranch="main",
    push_preview=true,
    forcepush=true,
)
