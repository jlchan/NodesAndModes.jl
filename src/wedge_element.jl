#####
##### 3D modes on pyramids
#####

"""
    basis(N, r, s, t)

Computes orthonormal basis of degree N at points (r,s,t)
"""
function basis(N,r,s,t,tol=1e-12)
    V_tri,Vr_tri,Vs_tri = Tri.basis(N,r,s)
    V_1D,Vt_1D = Line.basis(N,t)
    Np = (N+1)*(N+2)*(N+1)÷2
    V,Vr,Vs,Vt = ntuple(x->zeros(length(r),Np),4)
    id = 1
    for i = 1:size(V_tri,2)
        for j = 1:size(V_1D,2)
            V[:,id]  = @. V_tri[:,i]  * V_1D[:,j]
            Vr[:,id] = @. Vr_tri[:,i] * V_1D[:,j]
            Vs[:,id] = @. Vs_tri[:,i] * V_1D[:,j]
            Vt[:,id] = @. V_tri[:,i]  * Vt_1D[:,j]
            id += 1
        end
    end
    return V,Vr,Vs,Vt
end

"""
    equi_nodes(N)

Compute equispaced nodes of degree N.
"""
function equi_nodes(N)
    t1D = LinRange(-1,1,N+1)
    r_tri,s_tri = Tri.equi_nodes(N)
    r,t = vec.(meshgrid(r_tri,t1D))
    s,_ = vec.(meshgrid(s_tri,t1D))
    return r,s,t
end

"""
    nodes(N)

Computes interpolation nodes of degree N.
"""
function nodes(N)
    t1D,_ = gauss_lobatto_quad(0,0,N)
    r_tri,s_tri = Tri.nodes(N)
    r,t = vec.(meshgrid(r_tri,t1D))
    s,_ = vec.(meshgrid(s_tri,t1D))
    return r,s,t
end

"""
    quad_nodes(N)

Returns quadrature nodes and weights which exactly integrate degree 2N polynomials.
"""
function quad_nodes(N)
    r_tri,s_tri,w_tri = Tri.quad_nodes(2*N)
    t1D,wt1D = gauss_quad(0,0,N)
    wrs,wt = vec.(meshgrid(w_tri,wt1D))
    r,t = vec.(meshgrid(r_tri,t1D))
    s,_ = vec.(meshgrid(s_tri,t1D))
    return r,s,t,wrs.*wt
end
