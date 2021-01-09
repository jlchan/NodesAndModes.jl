#####
##### routines to compute warp-and-blend interpolation nodes for triangles, tets,
##### and pyramid elements. Wedges/hexes use tensor product nodes
#####

"""
    get_edge_list(elem::AbstractElemShape)

Returns list of edges for a specific element (elem = Tri(), Pyr(), or Tet()).
"""
get_edge_list(elem::AbstractElemShape)

get_edge_list(elem::Tri) = [1,2],[2,3],[3,1]
get_edge_list(elem::Tet) = [1,4],[4,3],[3,1],[1,2],[3,2],[4,2]
get_edge_list(elem::Pyr) = [1,2],[2,4],[3,4],[3,1],[1,5],[2,5],[3,5],[4,5]

get_vertices(elem::AbstractElemShape) = equi_nodes(elem,1)
get_vertex_fxns(elem::AbstractElemShape) =
    (rst...)->vandermonde(elem,1,rst...)/vandermonde(elem,1,equi_nodes(elem,1)...)

# assumes r1D in [-1,1], v1,v2 are vertices
function interp_1D_to_line(r1D,v1,v2)
    runit = @. (1+r1D)/2
    return map((a,b)->(b-a)*runit .+ a,v1,v2)
end

"""
    interp_1D_to_edges(elem::AbstractElemShape, r1D)

Interpolates points r1D to the edges of an element (elem = :Tri, :Pyr, or :Tet)
"""
function interp_1D_to_edges(elem::AbstractElemShape, r1D)
    v = get_vertices(elem)
    edges = get_edge_list(elem)
    edge_pts = ntuple(x->zeros(length(r1D),length(edges)),length(v))
    for (i,e) in enumerate(edges)
        edge_pts_i = interp_1D_to_line(r1D,getindex.(v,e[1]),getindex.(v,e[2]))
        setindex!.(edge_pts,edge_pts_i,:,i)
    end
    return vec.(edge_pts)
end

"""
    edge_basis(elem::AbstractElemShape, N, rst...)

returns generalized Vandermonde matrix evaluated using an edge basis.
"""
function edge_basis(elem::AbstractElemShape, N, rst...)
    vertices = get_vertices(elem)
    edges = get_edge_list(elem)
    vertex_functions = get_vertex_fxns(elem)
    return edge_basis(N, vertices, edges, (N,r)->basis(Line(),N,r), vertex_functions, rst...)
end

"""
    edge_basis(N, vertices, edges, basis, vertex_functions, rst...)

Computes edge basis given vertex functions and 1D basis.
"""
function edge_basis(N, vertices, edges, basis1D, vertex_functions, rst...)
    V1 = vertex_functions(rst...)
    if N < 2
        return V1
    end

    V = zeros(length(first(rst)),length(edges)*(N-1) + length(first(vertices)))
    V[:,1:size(V1,2)] .= V1 # initialize vertex functions

    id = size(V1,2)+1
    for e in edges
        # edge parametrization from paper:
        # "A Comparison of High-Order Interpolation Nodes for the Pyramid"
        r1D_edge = V1[:,e[1]]-V1[:,e[2]]
        V1D,_ = basis1D(N-2,r1D_edge)
        for i = 1:size(V1D,2)
            V[:,id] = V1D[:,i].*V1[:,e[1]].*V1[:,e[2]]
            id += 1
        end
    end
    return V
end

"""
    build_warped_nodes(elem::AbstractElemShape,N,r1D)

Computes degree N warp-and-blend interpolation nodes for elem = Tri(), Pyr(), or
Tet() based on the 1D node set "r1D". Returns a tuple "rst" containing arrays of
interpolation points.
"""
function build_warped_nodes(elem::AbstractElemShape,N,r1D)
    r1D_equi = equi_nodes(Line(),N)
    rst_edge_equi = interp_1D_to_edges(elem,r1D_equi)
    V_edge = edge_basis(elem,N,rst_edge_equi...)
    rst_edge = interp_1D_to_edges(elem,r1D)

    c = (x->V_edge\x).(rst_edge) # should be lsq solve

    rst_equi = equi_nodes(elem,N)
    return (x->edge_basis(elem,N,rst_equi...)*x).(c)
end


# function compute_lobatto_interp(elem::AbstractElemShape,N,r1D)
#     r1D_equi = equi_nodes(Line(),N)
#     rst_edge_equi = interp_1D_to_edges(elem,r1D_equi)
#     V_edge = edge_basis(elem,N,rst_edge_equi...)
#     rst_edge = interp_1D_to_edges(elem,r1D)
#
#     c = (x->V_edge\x).(rst_edge)
#     # @show ((x->V_edge*x).(c) .- rst_edge)
#     return c
# end
