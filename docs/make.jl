# using Pkg; pkg"activate .." # only need this for local builds?
using Documenter
using NodesAndModes

makedocs(
    sitename = "NodesAndModes.jl",
    repo = "https://github.com/jlchan/NodesAndModes.jl",
    modules=[NodesAndModes],
    pages = [
        "Home" => "index.md",
        "Line (1D) element" => "Line.md",
        "Triangular element" => "Tri.md",
        "Quadrilateral element" => "Quad.md",
        "Tetrahedral element" => "Tet.md",
        "Hexahedra element" => "Hex.md",
        "Wedge element" => "Wedge.md",
        "Pyramid element" => "Pyr.md",
        "Helper functions" => "nodes_and_modes_1D.md",
        "Authors" => "authors.md"
    ]
)

deploydocs(
    repo = "github.com/jlchan/NodesAndModes.jl",
)
