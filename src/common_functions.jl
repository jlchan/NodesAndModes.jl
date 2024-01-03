"""
    vandermonde(elem::AbstractElemShape, N, rst...)

Computes the generalized Vandermonde matrix V of degree N at points (r,s,t).
"""
vandermonde(elem::AbstractElemShape, N, rst...) = first(basis(elem,N,rst...))
vandermonde(elem::Line, N, r) = first(basis(elem,N,r)) # specialize for 1D

"""
    grad_vandermonde(elem::AbstractElemShape, N, rst...)

Computes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).
"""
grad_vandermonde(elem::AbstractElemShape, N, rst...) = basis(elem, N, rst...)[2:end]
grad_vandermonde(elem::Line, N, r) = last(basis(elem,N,r)) # specialize for 1D


"""
    nodes(elem::AbstractElemShape,N)

Computes interpolation nodes of degree N. Edge nodes coincide with (N+1)-point Lobatto points.
Default routine for elem = Tet(), Pyr(), Tri().

For Quad(), Hex(), Wedge() elements, nodes(...) returns interpolation points constructed
using a tensor product of lower-dimensional nodes. 
"""
nodes(elem::AbstractElemShape,N) = build_warped_nodes(elem, N, nodes(Line(),N))



"""
    basis(elem::AbstractElemShape, N, rst...)

Computes orthonormal basis of degree N at tuple of coordinate arrays (r,s,t).
"""
basis(elem::AbstractElemShape,N,rst...)


"""
    equi_nodes(elem::AbstractElemShape, N)

Compute equispaced nodes of degree N.
"""
equi_nodes(elem::AbstractElemShape, N)

"""
    quad_nodes(elem::AbstractElemShape, N)

Compute quadrature nodes and weights exact for (at least) degree 2N polynomials.
"""
quad_nodes(elem::AbstractElemShape, N)


"""
    stroud_quad_nodes(elem::AbstractElemShape,N)

Returns Stroud-type quadrature nodes and weights constructed from the tensor product
of (N+1)-point Gauss-Jacobi rules. Exact for degree 2N polynomials
"""
stroud_quad_nodes(elem::AbstractElemShape, N)
