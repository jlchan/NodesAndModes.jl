"""
    vandermonde(N,rst...)

Computes the generalized Vandermonde matrix V of degree N at points (r,s,t).
"""
vandermonde(N, rst...) = first(basis(N,rst...))

"""
    vandermonde(N,r)

Computes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case
"""
vandermonde(N, r) = first(basis(N,r)) # specialize for 1D

"""
    grad_vandermonde(N,rst...)

Computes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).
"""
grad_vandermonde(N, rst...) = basis(N,rst...)[2:end]

"""
    grad_vandermonde(N,r)

Computes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case
"""
grad_vandermonde(N, r) = last(basis(N,r)) # specialize for 1D
