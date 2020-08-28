module NodesAndModes

using LinearAlgebra
using SpecialFunctions

export greet
greet() = print("Hello World!")

include("nodes_and_modes.jl")

# export 1D routines by default (there's only one type of element in 1D)
export gauss_lobatto_quad, gauss_quad
export jacobiP, grad_jacobiP

# export 1D basis by default
export vandermonde_1D, grad_vandermonde_1D

"""
    vandermonde_1D(N, r)

Initialize the 1D Vandermonde matrix of order N Legendre polynomials at nodes r

# Examples
```jldoctest
"""
function vandermonde_1D(N, r)
    V1D = zeros(length(r), N+1)
    for j = 1:N+1
        V1D[:,j] = jacobiP(r[:], 0, 0, j-1)
    end
    return V1D
end

"""
    grad_vandermonde_1D(N, r)

Initialize the 1D Vandermonde matrix of order N Legendre polynomials at nodes r

# Examples
```jldoctest
"""
function grad_vandermonde_1D(N, r)
    V1D = zeros(length(r), N+1)
    for j = 1:N+1
        V1D[:,j] = grad_jacobiP(r[:], 0, 0, j-1)
    end
    return V1D
end




end # module
