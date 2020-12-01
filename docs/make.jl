using Pkg
pkg"activate .."

using Documenter
using NodesAndModes

makedocs(
    sitename = "NodesAndModes.jl",
    repo = "https://github.com/jlchan/NodesAndModes.jl",
    modules=[NodesAndModes],
    pages = [
        "Home" => "index.md",
        "Core 1D routines" => "nodes_and_modes_1D.md",
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
