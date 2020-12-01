"""
    vandermonde(N, r, s, t)

Computes generalized Vandermonde matrix.
"""
vandermonde(N, r, s, t) = first(basis(N,r,s,t))
"""
    grad_vandermonde(N, r, s, t)

Computes generalized Vandermonde-derivative matrices.
"""
grad_vandermonde(N, r, s, t) = basis(N,r,s,t)[2:4]

"""
    basis(N, r, s, t)

Computes orthonormal basis of degree N at points (r,s,t)

# Examples
```jldoctest
julia> N = 1; r,s,t = nodes(N);

julia> V,Vr,Vs,Vt = basis(N,r,s,t);

julia> V
 8×8 Array{Float64,2}:
 0.353553  -0.612372  -0.612372   1.06066  -0.612372   1.06066   1.06066  -1.83712
 0.353553  -0.612372   0.612372  -1.06066  -0.612372   1.06066  -1.06066   1.83712
 0.353553  -0.612372  -0.612372   1.06066   0.612372  -1.06066  -1.06066   1.83712
 ⋮                                                     ⋮
 0.353553   0.612372   0.612372   1.06066  -0.612372  -1.06066  -1.06066  -1.83712
 0.353553   0.612372  -0.612372  -1.06066   0.612372   1.06066  -1.06066  -1.83712
 0.353553   0.612372   0.612372   1.06066   0.612372   1.06066   1.06066   1.83712
```
"""
function basis(N, r, s, t)
    Np = convert(Int,(N+1)^3)
    sk = 1
    V,Vr,Vs,Vt = ntuple(x->zeros(length(r), Np),4)
    for i=0:N
        for j=0:N
            for k=0:N
                V[:,sk]  = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t, 0, 0, k)
                Vr[:,sk] = grad_jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vs[:,sk] = jacobiP(r, 0, 0, i).*grad_jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vt[:,sk] = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*grad_jacobiP(t,0,0,k)
                sk += 1
            end
        end
    end

    return V, Vr, Vs, Vt
end

# ===================================================

"""
    nodes(N)

Computes optimized interpolation nodes of degree N.

# Examples
```jldoctest
julia> N = 1; r,s,t = nodes(N);

julia> [r s t]
 8×3 Array{Float64,2}:
 -1.0  -1.0  -1.0
 -1.0   1.0  -1.0
  1.0  -1.0  -1.0
  1.0   1.0  -1.0
 -1.0  -1.0   1.0
 -1.0   1.0   1.0
  1.0  -1.0   1.0
  1.0   1.0   1.0
```
"""
function nodes(N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    equi_nodes(N)

Compute equispaced nodes of degree N.

# Examples
```jldoctest
julia> r,s,t = equi_nodes(1);
```
"""
function equi_nodes(N)
    r1D = LinRange(-1,1,N+1)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    quad_nodes(N)

Compute quadrature nodes and weights which exactly integrate degree 2N polynomials.

# Examples
```jldoctest
julia> r,s,t,w = quad_nodes(2)
```
"""
function quad_nodes(N)
    r1D,w1D = gauss_quad(0,0,N)
    r,s,t = vec.(meshgrid(r1D,r1D,r1D))
    wr,ws,wt = vec.(meshgrid(w1D,w1D,w1D))
    w = @. wr*ws*wt
    return r,s,t,w
end
