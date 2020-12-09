"""
    vandermonde(elem::ElementShape, N, rst...)

Computes the generalized Vandermonde matrix V of degree N at points (r,s,t).
"""
vandermonde(elem::ElementShape, N, rst...) = first(basis(elem,N,rst...))
vandermonde(elem::ElementShape, N, r) = first(basis(elem,N,r)) # specialize for 1D

"""
    grad_vandermonde(elem::ElementShape, N, rst...)

Computes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).
"""
grad_vandermonde(elem::ElementShape, N, rst...) = basis(elem, N, rst...)[2:end]
grad_vandermonde(elem::ElementShape, N, r) = last(basis(elem,N,r)) # specialize for 1D


"""
    nodes(elem::ElementShape,N)

Computes interpolation nodes of degree N. Edge nodes coincide with (N+1)-point Lobatto points.
Default routine for elem = Tet(), Pyr(), Tri().

For Quad(), Hex(), Wedge() elements, nodes(...) returns interpolation points constructed
using a tensor product of lower-dimensional nodes. 
"""
nodes(elem::ElementShape,N) = build_warped_nodes(elem,N,nodes(Line(),N))



"""
    basis(elem::ElementShape, N, rst...)

Computes orthonormal basis of degree N at tuple of coordinate arrays (r,s,t).
"""
basis(elem::ElementShape,N,rst...)


"""
    equi_nodes(elem::ElementShape, N)

Compute equispaced nodes of degree N.
"""
equi_nodes(elem::ElementShape, N)

"""
    quad_nodes(elem::ElementShape, N)

Compute quadrature nodes and weights exact for (at least) degree 2N polynomials.
"""
quad_nodes(elem::ElementShape, N)


"""
    stroud_quad_nodes(elem::ElementShape,N)

Returns Stroud-type quadrature nodes and weights constructed from the tensor product
of (N+1)-point Gauss-Jacobi rules. Exact for degree 2N polynomials
"""
stroud_quad_nodes(elem::ElementShape,N)
