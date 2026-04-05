using Documenter
using H3X

makedocs(;
    sitename = "H3X.jl",
    modules = [H3X],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://charleskawczynski.github.io/H3X.jl",
    ),
    pages = [
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
    warnonly = true,
)

deploydocs(;
    repo = "github.com/charleskawczynski/H3X.jl.git",
    devbranch = "main",
    push_preview = true,
)
