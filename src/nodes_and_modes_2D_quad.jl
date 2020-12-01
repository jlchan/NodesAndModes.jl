"""
    basis(N,r,s)

Computes orthonormal basis of degree N at coordinates (r,s)

# Examples
```jldoctest
julia> N = 1; r,s = nodes(N);

julia> V,Vr,Vs = basis(N,r,s);

julia> V
 4×4 Array{Float64,2}:
 0.5  -0.866025  -0.866025   1.5
 0.5  -0.866025   0.866025  -1.5
 0.5   0.866025  -0.866025  -1.5
 0.5   0.866025   0.866025   1.5
```
"""
function basis(N,r,s)
    Np = convert(Int,(N+1)*(N+1))
    sk = 1
    V,Vr,Vs = ntuple(x->zeros(length(r), Np),3)
    for i=0:N
        for j=0:N
            V[:,sk]  = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j)
            Vr[:,sk] = grad_jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j)
            Vs[:,sk] = jacobiP(r, 0, 0, i).*grad_jacobiP(s, 0, 0, j)
            sk += 1
        end
    end
    return V,Vr,Vs
end

# ===================================================

"""
    nodes(N)

Compute optimized interpolation nodes of degree N

# Examples
```jldoctest
julia> N = 1; r,s = nodes(N);

julia> r
 4-element Array{Float64,1}:
 -1.0
  1.0
 -1.0
  1.0

julia> s
 4-element Array{Float64,1}:
 -1.0
 -1.0
  1.0
  1.0
```
"""

function nodes(N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
    s,r = meshgrid(r1D)
    return r[:], s[:]
end

"""
    equi_nodes(N)

Compute equispaced nodes of degree N.

# Examples
```jldoctest
julia> N = 1; r,s = equi_nodes(N);
```
"""

function equi_nodes(N)
    r1D = LinRange(-1,1,N+1)
    s,r = meshgrid(r1D)
    return r[:], s[:]
end

"""
    quad_nodes(N)

Compute quadrature nodes and weights of degree N

# Examples
```jldoctest
julia> N = 1; r,s,w = quad_nodes(N);

julia> [r s w]
 4×3 Array{Float64,2}:
 -0.57735  -0.57735  1.0
  0.57735  -0.57735  1.0
 -0.57735   0.57735  1.0
  0.57735   0.57735  1.0
```
"""
function quad_nodes(N)
    r1D,w1D = gauss_quad(0,0,N)
    s,r = meshgrid(r1D)
    ws,wr = meshgrid(w1D)
    w = @. wr*ws
    return r[:], s[:], w[:]
end
