#####
##### 3D modes on pyramids
#####

vandermonde(N, r, s, t) = first(basis(N,r,s,t))
grad_vandermonde(N, r, s, t) = basis(N,r,s,t)[2:4]

"""
    basis(N,r,s,t,tol=1e-12)

Computes orthonormal semi-nodal basis on the biunit pyramid element.

Warning: nodal derivative matrices may contain errors for nodes at t = 1.
A way to avoid this is to use weak differentiation matrices computed using
quadrature rules with only interior nodes.

# Examples
```jldoctest
julia> N = 1; r,s,t,w = quad_nodes(N);

julia> V,Vr,Vs,Vt = basis(N,r,s,t);

julia> V
 8×5 Array{Float64,2}:
 -0.237249  1.38743   0.0       0.0       0.0
  1.05375   0.720759  0.0       0.0       0.0
 -0.237249  0.0       0.0       1.38743   0.0
  1.05375   0.0       0.0       0.720759  0.0
 -0.237249  0.0       1.38743   0.0       0.0
  1.05375   0.0       0.720759  0.0       0.0
 -0.237249  0.0       0.0       0.0       1.38743
  1.05375   0.0       0.0       0.0       0.720759
```
"""
function basis(N,r,s,t,tol=1e-12)

    # convert to abc
    a = @. 2*(r+1)/(1-t)-1
    b = @. 2*(s+1)/(1-t)-1
    c = t

    # fix top point
    ids = @. abs(1-t) < tol
    t[ids] .= 1-tol # WARNING: hacky (note factor of (1-c) = (1-t) in "pc" below)
    a[ids] .= -1
    b[ids] .= -1

    # change of vars from a to b
    dadr = @. 2/(1-t)
    dbds = @. 2/(1-t)
    dadt = @. (2*r + 2)/(1 - t)^2
    dbdt = @. (2*s + 2)/(1 - t)^2
    dcdt = 1

    Np = (N+1)*(N+2)*(2*N+3) ÷ 6
    V,Vr,Vs,Vt = ntuple(x->zeros(length(r),Np),4)

    ind = 1
    for k = 0:N
        # make nodal quad basis
        bq,aq,wab = Quad.quad_nodes(k)
        VDM,_ = Quad.basis(k,aq,bq)
        VDMab,Va,Vb = Quad.basis(k,a,b)
        Vab,DVa,DVb = map(A->A/VDM,(VDMab,Va,Vb))

        CNk = (N+2) / (2^(2*k+2)*(2*k+3))
        pNk = jacobiP(c,2*k+3,0,N-k)
        pc = @. ((1-c)/2)^k*pNk
        dpNk = grad_jacobiP(c,2*k+3,0,N-k)
        dpc = @. ((1-c)/2)^k*dpNk - (k/2^k)*((1 - c)^(k - 1))*pNk

        for ij = 1:(k+1)^2
            scale = 1/sqrt(CNk*wab[ij]) # normalization
            V[:,ind]  = @. scale*Vab[:,ij]*pc
            Vr[:,ind] = @. scale*(DVa[:,ij]*dadr)*pc
            Vs[:,ind] = @. scale*(DVb[:,ij]*dbds)*pc
            Vt[:,ind] = @. scale*((DVa[:,ij]*dadt + DVb[:,ij]*dbdt)*pc + Vab[:,ij]*dpc*dcdt)
            ind += 1
        end
    end
    return V,Vr,Vs,Vt
end

"""
    abctorst(a,b,c)

Converts from Stroud coordinates (a,b,c) on [-1,1]^3 to reference element
coordinates (r,s,t).
"""
function abctorst(a,b,c)
    r = @. .5*(1+a)*(1-c) - 1
    s = @. .5*(1+b)*(1-c) - 1
    t = @. c
    return r,s,t
end

"""
    nodes(N)

Computes interpolation nodes of degree N.
"""
function nodes(N)
    r1D,_ = gauss_lobatto_quad(0,0,N)
    return build_warped_nodes(N,:Pyr,r1D)
end

"""
    equi_nodes(N)

Computes equispaced nodes of degree N.
"""
function equi_nodes(N)
    Np = (N+1)*(N+2)*(2*N+3)/6
    a,b,c = ntuple(x->Float64[],3)
    c1D = LinRange(-1,1,N+1)
    for k = N:-1:0
        if k==0
            r1D = [.5]
        else
            r1D = LinRange(-1,1,k+1)
        end
        ak,bk = vec.(meshgrid(r1D,r1D))
        ck = c1D[N+1-k]*one.(ak)
        append!.((a,b,c),(ak,bk,ck))
    end
    return abctorst(a,b,c)
end

"""
    quad_nodes(N)

Computes quadrature nodes and weights which are exact for degree N polynomials.
"""
function quad_nodes(N)
    return stroud_quad_nodes(N)
end

"""
    stroud_quad_nodes(N)

Computes quadrature nodes and weights which are exact for degree N polynomials.
"""
function stroud_quad_nodes(N)

    a1D, w1D = gauss_quad(0,0,N)
    c1D, wc1D = gauss_quad(2,0,N)

    a,c,b = vec.(meshgrid(a1D, c1D, a1D))
    wa,wc,wb = vec.(meshgrid(w1D,wc1D,w1D))
    w = @. wa*wb*wc
    w = (8/3)*w./sum(w); # scale by b*h/3 = volume of pyr. b = 4, h = 2 for biunit

    return abctorst(a,b,c)...,w
end
