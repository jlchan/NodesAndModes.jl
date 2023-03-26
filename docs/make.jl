using Documenter
using NodesAndModes

makedocs(
    sitename = "NodesAndModes.jl",
    repo = "https://github.com/jlchan/NodesAndModes.jl",
    modules=[NodesAndModes],
    pages = [
        "Home" => "index.md",
        "Index" => "function_index.md",
        "Authors" => "authors.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/NodesAndModes.jl",
)
