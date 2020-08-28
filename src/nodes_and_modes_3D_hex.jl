"""
    vandermonde_3D(N, r)

Initialize the 3D Vandermonde matrix of order N "Legendre" polynomials at
nodes (r,s,t)

# Examples
N = 2
r,s,t = nodes_3D(N)
V = vandermonde_3D(N, r, s, t)

```jldoctest
"""
function vandermonde_3D(N, r, s, t)

    Np = convert(Int,(N+1)^3)
    sk = 1
    V = zeros(length(r), Np);
    for i=0:N
        for j=0:N
            for k=0:N
                V[:,sk] = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t, 0, 0, k)
                sk += 1
            end
        end
    end
    return V
end

"""
    grad_vandermonde_3D(N, r, s, t)

# Examples
N = 2
r,s,t = nodes_3D(N)
V = grad_vandermonde_3D(N, r, s, t)
```jldoctest
"""
function grad_vandermonde_3D(N, r, s, t)

    Np = convert(Int,(N+1)^3)
    sk = 1
    Vr = zeros(length(r), Np);
    Vs = zeros(length(r), Np);
    Vt = zeros(length(r), Np);
    for i=0:N
        for j=0:N
            for k=0:N
                Vr[:,sk] = grad_jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vs[:,sk] = jacobiP(r, 0, 0, i).*grad_jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vt[:,sk] = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*grad_jacobiP(t,0,0,k)
                sk += 1
            end
        end
    end

    return Vr, Vs, Vt
end

# ===================================================

"""
    nodes_3D(N)

# Examples
r,s,t = nodes_3D(N)
```jldoctest
"""
function nodes_3D(N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    equi_nodes_3D(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
r,s,t = equi_nodes_3D(N)
```jldoctest
"""

function equi_nodes_3D(N)
    r1D = LinRange(-1,1,N+1)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    quad_nodes_3D(N)
# Examples
N = 2
rq,sq,tq,wq = quad_nodes_3D(N)
```jldoctest
"""
function quad_nodes_3D(N)
    r1D,w1D = gauss_quad(0,0,N)
    r,s,t = vec.(meshgrid(r1D,r1D,r1D))
    wr,ws,wt = vec.(meshgrid(w1D,w1D,w1D))
    w = @. wr*ws*wt
    return r,s,t,w
end
