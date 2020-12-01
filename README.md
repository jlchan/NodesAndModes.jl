# NodesAndModes
[![Build Status](https://travis-ci.com/jlchan/NodesAndModes.jl.svg?branch=master)](https://travis-ci.com/jlchan/NodesAndModes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/NodesAndModes.jl?svg=true)](https://ci.appveyor.com/project/jlchan/NodesAndModes-jl)
[![Build status](https://github.com/jlchan/NodesAndModes.jl/workflows/CI/badge.svg)](https://github.com/jlchan/NodesAndModes.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/NodesAndModes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/NodesAndModes.jl)

Package to compute nodes and modes (orthogonal polynomials) for nodal DG methods on various reference elements. Codes inspired by the Matlab codes for the book "Nodal Discontinuous Galerkin methods" by Hesthaven and Warburton (2007) and high order interpolation nodes on triangles, tetrahedra, and pyramids using the "Interpolatory Warp and Blend" procedure from [Chan and Warburton 2015](https://epubs.siam.org/doi/abs/10.1137/141000105).

Original port from Matlab to Julia by [Yimin Lin](https://github.com/yiminllin).

This package is registered and can be installed via `] add NodesAndModes`. Julia v1.4 and above required.

## Example usage

To compute a 1D Vandermonde matrix using Gauss-Lobatto points and orthonormal polynomials.
```
using NodesAndModes.Line

N = 2
r,w = gauss_lobatto_quad(0,0,N)
V = vandermonde(N,r)
```

To compute a 2D Vandermonde matrix from Warp & Blend points (see [Warburton 2006](http://dx.doi.org/10.1007/s10665-006-9086-6)) and orthonormal polynomials on the triangle (with coordinates `r,s`)
```
using NodesAndModes.Tri

N = 2
r,s = nodes(N)
V = vandermonde(N,r,s)
```
Differentiation matrices `Dr` and `Ds` can be computed via
```
using NodesAndModes.Tri

N = 2
r,s = nodes(N)
V = vandermonde(N,r,s)
Dr,Ds = (A->A/V).(grad_vandermonde(N,r,s))
```

Nodes and modes (and mode derivatives) are available for 1D elements, 2D Tri and Quad elements, and 3D Hex, Tet, Wedge, and Pyr elements.

Each module exports the following functions:
- `vandermonde`: computes Vandermonde matrix constructed using orthogonal polynomials on the reference element.
- `grad_vandermonde`: computes matrix whose columns are derivatives of orthogonal polynomials on the reference element.
- `nodes`, `equi_nodes`: computes (non-uniform) interpolation nodes and equispaced nodes on the reference element.
- `quad_nodes`: computes quadrature nodes and weights on the reference element. 
