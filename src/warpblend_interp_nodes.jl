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
        return [1,4],[4,3],[3,1],[1,2],[3,2],[4,2]
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
