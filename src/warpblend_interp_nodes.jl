#####
##### routines to compute warp-and-blend interpolation nodes for triangles, tets,
##### and pyramid elements. Wedges/hexes use tensor product nodes
#####

"function get_edge_list(elem::Symbol)
    elem = :Tri, :Pyr, or :Tet.

    Example: edges = get_vertices_edges(:Tri)
"
function get_edge_list(elem::Symbol)
    if elem==:Tri
        return [1,2],[2,3],[3,1]
    elseif elem==:Pyr
        return [1,2],[2,4],[3,4],[3,1],[1,5],[2,5],[3,5],[4,5]
    elseif elem==:Tet
        return nothing # yet...
    end
end

# dispatch symbols to module
get_vertices(elem::Symbol) = get_vertices(@eval $elem)
get_vertex_fxns(elem::Symbol) = get_vertex_fxns(@eval $elem)
get_equi_nodes(N,elem::Symbol) = get_equi_nodes(N,@eval $elem)

# call actual module
get_vertices(elem::Module) = elem.equi_nodes(1)
get_vertex_fxns(elem::Module) =
    (rst...)->elem.vandermonde(1,rst...)/elem.vandermonde(1,elem.equi_nodes(1)...)
get_equi_nodes(N,elem::Module) = elem.equi_nodes(N)

# assumes r1D in [-1,1], v1,v2 are vertices
function interp_1D_to_line(r1D,v1,v2)
    runit = @. (1+r1D)/2
    return map((a,b)->(b-a)*runit .+ a,v1,v2)
end

function interp_1D_to_edges(r1D,elem::Symbol)
    v = get_vertices(elem)
    edges = get_edge_list(elem)
    edge_pts = ntuple(x->zeros(length(r1D),length(edges)),length(v))
    for (i,e) in enumerate(edges)
        edge_pts_i = interp_1D_to_line(r1D,getindex.(v,e[1]),getindex.(v,e[2]))
        setindex!.(edge_pts,edge_pts_i,:,i)
    end
    return vec.(edge_pts)
end

"
function edge_basis(N, elem::Symbol, rst...)
    build VDM matrix for edge basis
"
function edge_basis(N, elem::Symbol, rst...)
    vertices = get_vertices(elem)
    edges = get_edge_list(elem)
    vertex_functions = get_vertex_fxns(elem)
    return edge_basis(N, vertices, edges, Line.basis, vertex_functions, rst...)
end

function edge_basis(N, vertices, edges, basis, vertex_functions, rst...)
    V1 = vertex_functions(rst...)
    if N<2
        return V1
    end
    V = zeros(length(first(rst)),length(edges)*(N-1) + length(first(vertices)))
    V[:,1:size(V1,2)] .= V1
    id = size(V1,2)+1
    for e in edges
        # unit vector along which the edge runs
        dv = map(x->x[e[2]]-x[e[1]],vertices)
        dv = dv ./ sqrt(sum(dv.^2)) # normalize

        # eval 1D edge basis with edge parametrization
        r1D_edge = sum(rst .* dv)
        V1D,_ = Line.basis(N-2,r1D_edge)
        for i = 1:size(V1D,2)
            V[:,id] = V1D[:,i].*V1[:,e[1]].*V1[:,e[2]]
            id += 1
        end
    end
    return V
end

"
function build_warped_nodes(N,elem::Symbol,r1D)

Input:
    N:      polynomial degree
    elem:   :Tri, :Pyr, or :Tet.
    r1D:    1D node set for warping

Output:
    rst:    tuple containing arrays of interpolation points
"
function build_warped_nodes(N,elem::Symbol,r1D)
    r1D_equi = collect(LinRange(-1,1,N+1))
    rst_edge_equi = interp_1D_to_edges(r1D_equi,elem)
    V_edge = edge_basis(N,elem,rst_edge_equi...)
    rst_edge = interp_1D_to_edges(r1D,elem)

    c = (x->V_edge\x).(rst_edge)

    rst_equi = get_equi_nodes(N,elem)
    return (x->edge_basis(N,elem,rst_equi...)*x).(c)
end













"""
    xytors(x, y)

Transfer from equilateral triangle coordinate (a,b) to reference element
coordinate (r,s)

# Examples
```jldoctest
"""
function xytors(x, y)
    L₁ = @. (sqrt(3.0)y+1.0)/3.0
    L₂ = @. (-3.0x-sqrt(3.0)y+2.0)/6.0
    L₃ = @. (3.0x-sqrt(3.0)y+2.0)/6.0
    r = -L₂+L₃-L₁
    s = -L₂-L₃+L₁
    return r, s
end

"""
    function warp_factor(N, rout)

Compute warp factor for triangular interp nodes.

# Examples
```jldoctest
"""

function warp_factor(N, rout)
    LGLr,w1D = gauss_lobatto_quad(0, 0, N)
    # LGLr = transpose(LGLr[:])
    req = LinRange(-1,1,N+1)
    Veq = Line.vandermonde(N, req)
    Nr = length(rout)
    Pmat = zeros(N+1, Nr)
    for i = 1:N+1
        Pmat[i,:] = transpose(jacobiP(rout, 0, 0, i-1))
    end
    Lmat = transpose(Veq)\Pmat
    # warp = transpose(Lmat)*transpose(LGLr.-req)
    warp = transpose(Lmat)*(LGLr.-req[:])
    zerof = @. (abs(rout) < 1.0 - 1.0e-10)
    sf = @. 1.0 - (zerof*rout)^2
    warp = @. warp/sf + warp*(zerof-1)
    return warp
end


"""
    wb_nodes_tri(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function wb_nodes_tri(N)

    Np = convert(Int,(N+1)*(N+2)/2)

    alpopt = [0.0 0.0 1.4152 0.1001 0.2751 0.98 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.577 1.6223 1.6258]
    α = N < 16 ? alpopt[N] : 5/3

    L1 = zeros(Np, 1)
    L2 = zeros(Np, 1)
    L3 = zeros(Np, 1)
    sk = 1
    for n = 1:N+1
        for m = 1:N+2-n
            L1[sk] = (n-1)/N
            L3[sk] = (m-1)/N
            sk += 1;
        end
    end
    L2 = @. 1.0-L1-L3
    x = -L2+L3
    y = (-L2-L3+2*L1)/sqrt(3.0)
    blend1 = 4.0*L2.*L3
    blend2 = 4.0*L1.*L3
    blend3 = 4.0*L1.*L2
    warpf1 = warp_factor(N, L3-L2)
    warpf2 = warp_factor(N, L1-L3)
    warpf3 = warp_factor(N, L2-L1)
    warp1 = @. blend1*warpf1*(1+(α*L1).^2)
    warp2 = @. blend2*warpf2*(1+(α*L2).^2)
    warp3 = @. blend3*warpf3*(1+(α*L3).^2)

    x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3
    y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3

    r, s = xytors(x, y)
    return r[:], s[:]
end
