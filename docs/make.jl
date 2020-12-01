using Documenter, NodesAndModes

Documenter.makedocs(
    sitename = "jl documentation",
    repo = "https://github.com/jlchan/jl",
    format = Documenter.HTML(),
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
