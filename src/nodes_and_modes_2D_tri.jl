#####
##### 2D modes on triangles
#####

"""
    simplex_2D(a, b, i, j)

Evaluate 2D "Legendre" basis phi_ij at (a,b) returned by rstoab

# Examples
```jldoctest
"""
function simplex_2D(a, b, i, j)
    h₁ = jacobiP(a, 0, 0, i)
    h₂ = jacobiP(b, 2i+1, 0, j)
    return @. sqrt(2.0)*h₁*h₂*(1-b)^i
end

"""
    grad_simplex_2D(a, b, id, jd)

Evalute the partial derivative w.r.t. (r,s) of 2D "Legendre" polynomial with
index, or order, (id,jd) at (a,b)

# Examples
```jldoctest
"""

function grad_simplex_2D(a, b, id, jd)
    fa = jacobiP(a, 0, 0, id)
    gb = jacobiP(b, 2*id+1, 0, jd)
    dfa = grad_jacobiP(a, 0, 0, id)
    dgb = grad_jacobiP(b, 2*id+1, 0, jd)

    dmodedr = @. dfa*gb
    if id > 0
        dmodedr = dmodedr.*((0.5*(1 .- b)).^(id-1))
    end

    dmodeds = @. dfa*(gb*(0.5*(1+a)))
    if id > 0
        dmodeds = dmodeds.*((0.5*(1 .- b)).^(id-1))
    end

    tmp = @. dgb*((0.5*(1-b))^id)
    if id > 0
        tmp = tmp-0.5*id*gb.*((0.5*(1 .- b)).^(id-1))
    end
    dmodeds = dmodeds+fa.*tmp

    dmodedr = 2^(id+0.5)*dmodedr
    dmodeds = 2^(id+0.5)*dmodeds
    return dmodedr, dmodeds
end


"""
    rstoab(r, s)

Transfer from reference element coordinate (r,s) to polynomial basis evaluation
coordinate (a,b)

# Examples
```jldoctest
"""
function rstoab(r, s)
    a = zeros(length(r), 1)
    for n = 1:length(r)
        if s[n] != 1
            a[n] = 2*(1+r[n])/(1-s[n])-1
        else
            a[n] = -1
        end
    end
    return a, s
end


"""
    equi_nodes_tri(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function equi_nodes_tri(N)
    Np = convert(Int,(N+1)*(N+2)/2)

    r = zeros(Np,1)
    s = zeros(Np,1)

    r1D = LinRange(-1,1,N+1)
    sk = 1
    for i = 0:N
        for j = 0:N-i
            r[sk] = r1D[i+1]
            s[sk] = r1D[j+1]
            sk += 1
        end
    end

    return r[:], s[:]
end


"""
    quad_nodes_tri(N)

Returns quadrature nodes which exactly integrate degree N polynomials

# Examples
```jldoctest
"""

function quad_nodes_tri(N)

    if N < 28
        rsw = readdlm(string(@__DIR__,"/QuadratureData/quad_nodes_tri_N", N, ".txt"),' ', Float64, '\n')
        r = rsw[:,1]
        s = rsw[:,2]
        w = rsw[:,3]
    else
        cubN = convert(Int,ceil((N+1)/2))
        r,s,w = stroud_quad_nodes(cubN)
    end

    return r[:], s[:], w[:]
end

function stroud_quad_nodes(N)
    cubA,cubWA = gauss_quad(0,0,N)
    cubB,cubWB = gauss_quad(1,0,N)

    cubA = ones(N+1,1)*cubA'
    cubB = cubB*ones(1,N+1)

    r = @. 0.5*(1+cubA)*(1-cubB)-1
    s = cubB
    w = 0.5*cubWB*(cubWA')
    return vec.((r,s,w))
end


"""
    nodes(N)

# Examples
```jldoctest
"""
function nodes(N)
    return wb_nodes_tri(N)
end

"""
    equi_nodes(N)

# Examples
```jldoctest
"""
function equi_nodes(N)
    equi_nodes_tri(N)
end


"""
    quad_nodes(N)
    returns quadrature nodes which exactly integrate degree 2N polynomials

# Examples
```jldoctest
"""

function quad_nodes(N)
    r,s,w = quad_nodes_tri(2*N)
    return r,s,w
end

"""
function basis(N,r,s)
    basis: computes Vandermonde matrices for a basis + its derivatives
"""

function basis(N,r,s)
    Np = convert(Int,(N+1)*(N+2)/2)
    V2D,V2Dr,V2Ds = ntuple(x->zeros(length(r), Np),3)
    a, b = rstoab(r, s)
    sk = 1
    for i = 0:N
        for j = 0:N-i
            V2D[:,sk] = simplex_2D(a,b,i,j)
            V2Dr[:,sk], V2Ds[:,sk] = grad_simplex_2D(a, b, i, j)
            sk = sk+1
        end
    end
    return V2D,V2Dr,V2Ds
end

vandermonde(N, r, s) = first(basis(N,r,s))
grad_vandermonde(N, r, s) = basis(N,r,s)[2:3]
