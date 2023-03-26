# NodesAndModes

`NodesAndModes.jl` is a package to compute nodes (interpolation and quadrature points) and modes (orthogonal polynomials) on various reference elements for use in high order finite element and nodal discontinuous Galerkin (DG) methods.

The codes inspired by the Matlab codes for the book "Nodal Discontinuous Galerkin methods" by Hesthaven and Warburton (2007) and high order interpolation nodes on triangles, tetrahedra, and pyramids using the "Interpolatory Warp and Blend" procedure from [Chan and Warburton 2015](https://epubs.siam.org/doi/abs/10.1137/141000105).

## Installation

This package is registered and can be installed via `] add NodesAndModes`. Julia v1.3 is required.

## Package organization

`NodesAndModes.jl` supports seven types of elements in 1D, 2D, and 3D.
- Line (1D lines/intervals)
- Tri (2D triangles)
- Quad (2D quadrilaterals)
- Tet (3D tetrahedra)
- Hex (3D hexahedra)
- Wedge (3D wedges/prisms)
- Pyr (3D pyramids)

Each module exports the following functions:
- `basis`: computes Vandermonde matrix (columns are evaluations of orthogonal polynomials at different points) and derivative matrices (columns are derivatives of orthogonal polynomials at different points) constructed using orthogonal polynomials on the reference element
- `nodes`: computes (non-uniform) interpolation nodes on the reference element. All interpolation nodes
- `quad_nodes`: computes quadrature nodes and weights on the reference element. `quad_nodes(N)` returns a quadrature rule exact for degree 2N polynomials (e.g., exact integration of the mass matrix).
- `equi_nodes`: computes equispaced nodes on the reference element (for plotting)

!!! note

    Each module also exports `vandermonde` and `grad_vandermonde` for similarity with the Matlab codes in Hesthaven/Warburton 2007. These just call the `basis` routine and return either Vandermonde or derivative matrices.

## Example usage

To compute a 1D Vandermonde matrix using Gauss-Lobatto points and orthonormal polynomials.
```julia
julia> using NodesAndModes
julia> N = 2
julia> r = nodes(Line(),N)
julia> V = vandermonde(Line(),N,r)
```

To compute a 2D triangular Vandermonde matrix from Warp & Blend points (see [Warburton 2006](http://dx.doi.org/10.1007/s10665-006-9086-6)) and orthonormal polynomials on the triangle (with coordinates `r,s`)
```julia
julia> using NodesAndModes
julia> N = 2
julia> r,s = nodes(Tri(),N)
julia> V = vandermonde(Tri(),N,r,s)
```
Nodal differentiation matrices `Dr` and `Ds` can be computed via
```julia
julia> using NodesAndModes
julia> N = 2
julia> r,s = nodes(Tri(),N)
julia> V,Vr,Vs = basis(Tri(),N,r,s)
julia> Dr,Ds = (A->A/V).((Vr,Vs))
```
such that `Dr*f(r,s) ≈ df/dr`.

A mass matrix `M` and weak differentation matrices `Qr,Qs` in finite element or DG methods can be computed using quadrature via
```julia
julia> using LinearAlgebra
julia> using NodesAndModes
julia> N = 2
julia> r,s = nodes(Tri(),N)
julia> V = vandermonde(Tri(),N,r,s)
julia> rq,sq,wq = quad_nodes(Tri(),N)
julia> Vq,Vrq,Vsq = (A->A/V).(basis(Tri(),N,rq,sq))
julia> M,Qr,Qs = (A->Vq'*diagm(wq)*A).((Vq,Vrq,Vsq))
```
