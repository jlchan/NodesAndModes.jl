# NodesAndModes
[![Build Status](https://travis-ci.com/jlchan/NodesAndModes.jl.svg?branch=master)](https://travis-ci.com/jlchan/NodesAndModes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/NodesAndModes.jl?svg=true)](https://ci.appveyor.com/project/jlchan/NodesAndModes-jl)
[![Codecov](https://codecov.io/gh/jlchan/NodesAndModes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/NodesAndModes.jl)

Package to compute nodes and modes (orthogonal polynomials) for nodal DG methods. Codes inspired by the Matlab codes for the book Nodal Discontinuous Galerkin methods by Hesthaven and Warburton (2007). 

## Example usage

To compute a Vandermonde matrix from Lobatto points and orthonormal polynomials.
```
using NodesAndModes

N = 2
r,w = gauss_lobatto_quad(0,0,N)
V = vandermonde_1D(N,r)
```

To compute a Vandermonde matrix from Warp & Blend points (see [Warburton 2006](http://dx.doi.org/10.1007/s10665-006-9086-6)) and orthonormal polynomials on the triangle
```
using NodesAndModes.Tri

N = 2
r,s = nodes_2D(N)
V = vandermonde_2D(N,r,s)
```

Nodes and modes (and mode derivatives) available for 1D elements, 2D Tri and Quad elements, and 3D Hex elements. 

## To-do list
- support for tetrahedra, prisms, and pyramids. 
