#####
##### 2D modes on triangles
#####

"""
    simplex_2D(a, b, i, j)

Evaluate 2D PKDO basis phi_ij at points (a,b) on the Duffy domain [-1,1]^2
"""
function simplex_2D(a, b, i, j)
    h₁ = jacobiP(a, 0, 0, i)
    h₂ = jacobiP(b, 2i+1, 0, j)
    return @. sqrt(2.0)*h₁*h₂*(1-b)^i # WARNING: (1-b)^i can blow up if N is too large (> 25)
end

"""
    grad_simplex_2D(a, b, id, jd)

Evalute the partial derivative w.r.t. (r,s) of 2D PKDO polynomials with
indices (id,jd) at (a,b)
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
    rstoab(r, s, tol = 1e-12)

Converts from reference bi-unit right triangle coordinate (r,s) to polynomial
basis evaluation coordinates (a,b) on the domain [-1,1]^2
"""
function rstoab(r, s, tol = 1e-12)
    a = zeros(length(r))
    for n in eachindex(r)
        if abs(s[n] - 1) > tol
            a[n] = 2 * (1 + r[n]) / (1 - s[n]) - 1
        else
            a[n] = -1
        end
    end
    b = copy(s)
    return a, b
end

rstoab(::Tri, r, s, kwargs...) = rstoab(r, s, kwargs...)

"""
    abtors(Tri(), r, s, tol = 1e-12)

Converts from polynomial basis evaluation coordinates (a,b) on the 
domain [-1,1]^2 to reference bi-unit right triangle coordinate (r,s).
"""
function abtors(::Tri, a, b)
    r = @. 0.5 * (a + 1) * (1 - b) - 1
    s = copy(b)
    return r, s
end

"""
    quad_nodes_tri(N)

Returns quadrature nodes (from Gimbutas and Xiao 2010) which exactly integrate degree N polynomials
"""
function quad_nodes_tri(N)

    if N < 51
        rsw::Matrix{Float64} = readdlm(string(@__DIR__, "/QuadratureData/quad_nodes_tri_N", N, ".txt"), ' ', Float64, '\n')
        r = rsw[:, 1]
        s = rsw[:, 2]
        w = rsw[:, 3]
    else
        cubN = ceil(Int, (N+1) / 2)
        r, s, w = stroud_quad_nodes(Tri(), cubN)
    end

    return r[:], s[:], w[:]
end

function stroud_quad_nodes(::Tri, N)
    cubA, cubWA = gauss_quad(0, 0, N)
    cubB, cubWB = gauss_quad(1, 0, N)

    cubA, cubB = vec.(meshgrid(cubA, cubB))

    r, s = abtors(Tri(), cubA, cubB)
    w = 0.5 * cubWB * (cubWA')
    return vec.((r, s, w))
end

function equi_nodes(::Tri, N)
    Np = (N+1) * (N+2) ÷ 2

    r = zeros(Np)
    s = zeros(Np)

    r1D = LinRange(-1, 1, N+1)
    sk = 1

    # order so that the nodes change fastest in the x-direction,
    # then the y-direction. For N=1, the nodes should be
    #   
    # 3
    # | \
    # 1---2

    for j = 0:N, i = 0:N-j 
        r[sk] = r1D[i+1]
        s[sk] = r1D[j+1]
        sk += 1
    end

    return r[:], s[:]
end


function quad_nodes(::Tri, N)    
    return quad_nodes_tri(2*N)
end

function basis(::Tri, N, r, s)
    Np = (N+1) * (N+2) ÷ 2
    V2D, V2Dr, V2Ds = ntuple(x -> zeros(length(r), Np), 3)
    a, b = rstoab(r, s)
    sk = 1
    for i = 0:N
        for j = 0:N-i
            V2D[:, sk] = simplex_2D(a, b, i, j)
            V2Dr[:, sk], V2Ds[:, sk] = grad_simplex_2D(a, b, i, j)
            sk = sk + 1
        end
    end
    return V2D, V2Dr, V2Ds
end
