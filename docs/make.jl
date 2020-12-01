using Documenter, NodesAndModes

Documenter.makedocs(
    sitename = "NodesAndModes.jl documentation",
    repo = "https://github.com/jlchan/NodesAndModes.jl",
    pages = [
    "Home" => "index.md",
    "Line submodule" => "Line.md",
    "Tri submodule" => "Tri.md",
    "Quad submodule" => "Quad.md",
    "Tet submodule" => "Tet.md",
    "Hex submodule" => "Hex.md",
    "Wedge submodule" => "Wedge.md",
    "Pyr submodule" => "Pyr.md",
    "Authors" => "authors.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/NodesAndModes.jl",
)
