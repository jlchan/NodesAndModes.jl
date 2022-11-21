#####
##### 3D modes on pyramids
#####

function basis(::Wedge, N, r, s, t)
    V_tri, Vr_tri, Vs_tri = basis(Tri(), N, r, s)
    V_1D, Vt_1D = basis(Line(), N, t)
    Np = (N+1) * (N+2) * (N+1) ÷ 2
    V, Vr, Vs, Vt = ntuple(x->zeros(length(r), Np), 4)
    id = 1
    for j in axes(V_1D, 2)
        for i in axes(V_tri, 2)
            @. V[:,id]  = V_tri[:,i]  * V_1D[:,j]
            @. Vr[:,id] = Vr_tri[:,i] * V_1D[:,j]
            @. Vs[:,id] = Vs_tri[:,i] * V_1D[:,j]
            @. Vt[:,id] = V_tri[:,i]  * Vt_1D[:,j]
            id += 1
        end
    end
    return V,Vr,Vs,Vt
end

function equi_nodes(::Wedge,N)
    t1D = LinRange(-1, 1, N+1)
    r_tri, s_tri = equi_nodes(Tri(), N)
    t, r = vec.(meshgrid(t1D, r_tri))
    _, s = vec.(meshgrid(t1D, s_tri))    
    return r,s,t
end


function nodes(::Wedge, N)
    t1D, _ = gauss_lobatto_quad(0, 0, N)
    r_tri, s_tri = nodes(Tri(), N)
    t, r = vec.(meshgrid(t1D, r_tri))
    _, s = vec.(meshgrid(t1D, s_tri))    
    return r,s,t
end


function quad_nodes(::Wedge, N)
    r_tri, s_tri, w_tri = quad_nodes(Tri(), N)
    t1D, wt1D = quad_nodes(Line(), N)
    wt, wrs = vec.(meshgrid(wt1D, w_tri))
    t, r = vec.(meshgrid(t1D, r_tri))
    _, s = vec.(meshgrid(t1D, s_tri))  
    return r, s, t, wrs .* wt
end
