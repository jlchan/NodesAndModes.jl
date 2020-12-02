# NodesAndModes
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlchan.github.io/NodesAndModes.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlchan.github.io/NodesAndModes.jl/dev)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/NodesAndModes.jl?svg=true)](https://ci.appveyor.com/project/jlchan/NodesAndModes-jl)
[![Build status](https://github.com/jlchan/NodesAndModes.jl/workflows/CI/badge.svg)](https://github.com/jlchan/NodesAndModes.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/NodesAndModes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/NodesAndModes.jl)

Package to compute nodes and modes (orthogonal polynomials) for nodal DG methods on various reference elements. Codes inspired by the Matlab codes for the book "Nodal Discontinuous Galerkin methods" by Hesthaven and Warburton (2007) and high order interpolation nodes on triangles, tetrahedra, and pyramids using the "Interpolatory Warp and Blend" procedure from [Chan and Warburton 2015](https://epubs.siam.org/doi/abs/10.1137/141000105).

Original port from Matlab to Julia by [Yimin Lin](https://github.com/yiminllin).

This package is registered and can be installed via `] add NodesAndModes`. Julia v1.4 and above required.

Nodes and modes (and mode derivatives) are available for 1D elements, 2D Tri and Quad elements, and 3D Hex, Tet, Wedge, and Pyr elements.

Each submodule (Line, Tri, Quad, Hex, Tet, Wedge, Pyr) exports the following functions:
- `vandermonde`: computes Vandermonde matrix constructed using orthogonal polynomials on the reference element.
- `grad_vandermonde`: computes matrix whose columns are derivatives of orthogonal polynomials on the reference element.
- `basis`: returns outputs from both `vandermonde` and `grad_vandermonde`
- `nodes`: computes (non-uniform) interpolation nodes on the reference element.
- `quad_nodes`: computes quadrature nodes and weights on the reference element.
- `equi_nodes`: computes equispaced nodes on the reference element (for plotting)
