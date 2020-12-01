"""
    vandermonde(N,r)

Computes the generalized Vandermonde matrix V of degree N at points r.
"""
vandermonde(N, rst...) = first(basis(N,rst...))
vandermonde(N, r) = first(basis(N,r)) # specialize for 1D

"""
    grad_vandermonde(N,r,s)

Computes the generalized Vandermonde matrix V of degree N at points r.
"""
grad_vandermonde(N, rst...) = basis(N,rst...)[2:end]
grad_vandermonde(N, r) = last(basis(N,r)) # specialize for 1D
