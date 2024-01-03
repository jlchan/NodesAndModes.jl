#####
##### 3D modes on pyramids
#####

"""
    basis(elem::Pyr,N,r,s,t,tol=1e-12)

Computes orthonormal semi-nodal basis on the biunit pyramid element.

Warning: nodal derivative matrices may contain errors for nodes at t = 1.
A way to avoid this is to use weak differentiation matrices computed using
quadrature rules with only interior nodes.
"""
function basis(elem::Pyr, N, r, s, t, tol = 1e3*eps())

    # convert to abc
    a = @. 2 * (r + 1) / (1 - t) - 1
    b = @. 2 * (s + 1) / (1 - t) - 1
    c = t

    # fix top point
    ids = @. abs(1 - t) < tol
    t[ids] .= 1 - tol # WARNING: hacky (note factor of (1-c) = (1-t) in "pc" below)
    a[ids] .= -1
    b[ids] .= -1

    # change of vars from a to b
    dadr = @. 2 / (1 - t)
    dbds = @. 2 / (1 - t)
    dadt = @. (2 * r + 2) / (1 - t)^2
    dbdt = @. (2 * s + 2) / (1 - t)^2
    dcdt = 1

    Np = (N + 1) * (N + 2) * (2 * N + 3) ÷ 6
    V, Vr, Vs, Vt = ntuple(x->zeros(length(r), Np), 4)

    ind = 1
    for k = 0:N
        # make nodal quad basis
        bq, aq, wab = quad_nodes(Quad(), k)
        VDM, _ = basis(Quad(), k, aq, bq)
        VDMab, Va, Vb = basis(Quad(), k, a, b)
        Vab, DVa, DVb = map(A->A / VDM, (VDMab, Va, Vb))

        CNk = (N+2) / (2^(2*k+2) * (2*k+3))
        pNk = jacobiP(c, 2*k+3, 0, N-k)
        pc = @. ((1-c) / 2)^k * pNk
        dpNk = grad_jacobiP(c, 2*k+3, 0, N-k)
        dpc = @. ((1-c)/2)^k * dpNk - (k/2^k) * ((1 - c)^(k - 1)) * pNk

        for ij = 1:(k+1)^2
            scale = 1 / sqrt(CNk * wab[ij]) # normalization
            @. V[:,ind]  = scale * Vab[:,ij] * pc
            @. Vr[:,ind] = scale * (DVa[:,ij] * dadr) * pc
            @. Vs[:,ind] = scale * (DVb[:,ij] * dbds) * pc
            @. Vt[:,ind] = scale * ((DVa[:,ij] * dadt + DVb[:,ij] * dbdt) * pc + Vab[:,ij] * dpc * dcdt)
            ind += 1
        end
    end
    return V, Vr, Vs, Vt
end

"""
    abctorst(elem::Pyr,a,b,c)

Converts from Stroud coordinates (a,b,c) on [-1,1]^3 to reference
element coordinates (r,s,t).
"""
function abctorst(elem::Pyr, a, b, c)
    r = @. .5 * (1 + a) * (1 - c) - 1
    s = @. .5 * (1 + b) * (1 - c) - 1
    t = @. c
    return r, s, t
end

"""
    rsttoabc(elem::Pyr,a,b,c)

Converts from reference element coordinates (r,s,t) to Stroud 
coordinates (a,b,c) on [-1,1]^3.
"""
function rsttoabc(::Pyr, r, s, t)
    a = @. 2 * (r + 1) / (1 - t) - 1 
    b = @. 2 * (s + 1) / (1 - t) - 1
    c = @. t
    return a, b, c
end


"""
    nodes(elem::Pyr,N)

Computes interpolation nodes of degree N. Edge nodes coincide with (N+1)-point Lobatto
points. Triangular face nodes coincide with Tri.nodes(N), quadrilateral face nodes
coincide with tensor product (N+1)-point Lobatto points.
"""
function nodes(elem::Pyr, N)

    if N == 1
        return equi_nodes(Pyr(), N)
    end

    # append quad face nodes
    r1D_equi = equi_nodes(Line(), N)
    rst_equi = interp_1D_to_edges(Pyr(), r1D_equi)
    quad_face_pts_equi = (vec.(meshgrid(r1D_equi[2:end-1]))..., -ones((N - 1) * (N - 1)))
    append!.(rst_equi, quad_face_pts_equi)

    # append quad face nodes
    r1D = nodes(Line(), N)
    rst_lobatto = interp_1D_to_edges(Pyr(), r1D)
    quad_face_pts = (vec.(meshgrid(r1D[2:end-1]))..., -ones((N - 1) * (N - 1)))
    append!.(rst_lobatto, quad_face_pts)

    # append quad face "bubble" basis functions
    function edge_face_pyr_basis(N, rst)
        V_edge = edge_basis(Pyr(), N, rst...)
        # WARNING: assumes first 4 vertices of Pyr = quad vertices
        # also assumes first 5 cols of V_edge = vertex functions for Pyr
        V_quad_face = zeros(length(first(rst)), (N - 1) * (N - 1))
        V_quad = vandermonde(Quad(), N-2, rst[1], rst[2])
        V_quad_bubble = @. V_edge[:,1] * V_edge[:,2] * V_edge[:,3] * V_edge[:,4]
        for i = 1:size(V_quad_face, 2)
            @. V_quad_face[:,i] = V_quad_bubble * V_quad[:,i]
        end
        return hcat(V_edge, V_quad_face)
    end

    # same interp problem as before
    V_edge_face = edge_face_pyr_basis(N, rst_equi)
    c = (x->V_edge_face \ x).(rst_lobatto) # should be lsq solve
    return (x->edge_face_pyr_basis(N, equi_nodes(elem, N)) * x).(c)
end

function equi_nodes(elem::Pyr, N)
    a, b, c = ntuple(x -> Float64[], 3)
    c1D = LinRange(-1, 1, N+1)
    for k = N:-1:0
        if k==0
            r1D = [.5]
        else
            r1D = LinRange(-1, 1, k + 1)
        end
        ak, bk = vec.(meshgrid(r1D, r1D))
        ck = c1D[N + 1 - k] * one.(ak)
        append!.((a, b, c), (ak, bk, ck))
    end
    return abctorst(elem, a, b, c)
end

quad_nodes(elem::Pyr, N) = stroud_quad_nodes(elem, N)

function stroud_quad_nodes(elem::Pyr, N)

    a1D, w1D = gauss_quad(0, 0, N)
    c1D, wc1D = gauss_quad(2, 0, N)

    a, c, b = vec.(meshgrid(a1D, c1D, a1D))
    wa, wc, wb = vec.(meshgrid(w1D, wc1D, w1D))
    w = @. wa * wb * wc
    w = (8 / 3) * w ./ sum(w) # scale by b*h/3 = volume of pyr. b = 4, h = 2 for biunit

    return abctorst(elem, a, b, c)..., w
end
