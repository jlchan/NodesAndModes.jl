"""
    vandermonde_2D(N, r)

Initialize the 2D Vandermonde matrix of order N "Legendre" polynomials at
nodes (r,s)

# Examples
```jldoctest
"""
function vandermonde_2D(N, r, s)

    Np = convert(Int,(N+1)*(N+1))
    sk = 1
    V = zeros(length(r), Np);
    for i=0:N
        for j=0:N
            V[:,sk] = jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j)
            sk += 1
        end
    end
    return V
end

"""
    gradV_2D(N, Np, r, s)

Initilize 2D gradient Vandermonde matrix of derivatives 2D "Legendre" polynomials
of order N at (r,s)

# Examples
```jldoctest
"""
function grad_vandermonde_2D(N, r, s)

    Np = convert(Int,(N+1)*(N+1))
    sk = 1
    V2Dr = zeros(length(r), Np);
    V2Ds = zeros(length(r), Np);
    for i=0:N
        for j=0:N
            V2Dr[:,sk] = grad_jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j)
            V2Ds[:,sk] = jacobiP(r, 0, 0, i).*grad_jacobiP(s, 0, 0, j)
            sk += 1
        end
    end

    return V2Dr, V2Ds
end

# ===================================================

"""
    nodes_2D(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function nodes_2D(N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
    s,r = meshgrid(r1D)
    return r[:], s[:]
end

"""
    equi_nodes_2D(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function equi_nodes_2D(N)
    r1D = LinRange(-1,1,N+1)
    s,r = meshgrid(r1D)

    return r[:], s[:]
end

"""
    quad_nodes_2D(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function quad_nodes_2D(N)
    r1D,w1D = gauss_quad(0,0,N)
    s,r = meshgrid(r1D)
    ws,wr = meshgrid(w1D)
    w = @. wr*ws
    return r[:], s[:], w[:]
end
