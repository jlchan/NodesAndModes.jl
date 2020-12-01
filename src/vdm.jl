"""
    vandermonde(N,r)

Computes the generalized Vandermonde matrix V of degree N at points r.

# Examples
```jldoctest
julia> N = 1; r = Line.nodes(N);  V = Line.vandermonde(N,r)
julia> N = 1; r,s = Tri.nodes(N); V = Tri.vandermonde(N,r,s)
```
"""
vandermonde(N, rst...) = first(basis(N,rst...))
vandermonde(N, r) = first(basis(N,r)) # specialize for 1D

"""
    grad_vandermonde(N,r,s)

Computes the generalized Vandermonde matrix V of degree N at points r.

# Examples
```jldoctest
julia> N = 1; r = Line.nodes(N);  Vr = Line.grad_vandermonde(N,r)
julia> N = 1; r,s = Tri.nodes(N); Vr = Tri.grad_vandermonde(N,r,s)
```
"""
grad_vandermonde(N, rst...) = basis(N,rst...)[2:end]
grad_vandermonde(N, r) = last(basis(N,r)) # specialize for 1D
