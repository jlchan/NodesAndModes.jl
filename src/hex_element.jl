"""
    basis(N, r, s, t)

Computes orthonormal basis of degree N at points (r,s,t).
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
"""
function nodes(N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    equi_nodes(N)

Compute equispaced nodes of degree N.
"""
function equi_nodes(N)
    r1D = LinRange(-1,1,N+1)
    return vec.(meshgrid(r1D,r1D,r1D))
end

"""
    quad_nodes(N)

Compute quadrature nodes and weights exact for degree 2N+1 polynomials. Uses a tensor product
Gauss quadrature rule.
"""
function quad_nodes(N)
    r1D,w1D = gauss_quad(0,0,N)
    r,s,t = vec.(meshgrid(r1D,r1D,r1D))
    wr,ws,wt = vec.(meshgrid(w1D,w1D,w1D))
    w = @. wr*ws*wt
    return r,s,t,w
end
