#####
##### routines to compute warp-and-blend interpolation nodes for triangles, tets,
##### and pyramid elements. Wedges/hexes use tensor product nodes
#####

"""
    get_edge_list(elem::AbstractElemShape)

Returns list of edges for a specific element (elem = Tri(), Pyr(), Hex(), or Tet()).
"""
get_edge_list(elem::AbstractElemShape)

get_edge_list(elem::Union{Tri, Quad}) = SVector{2}.(find_face_nodes(elem, equi_nodes(elem, 1)...))
get_edge_list(::Tet)  = SVector(1, 4), SVector(4, 3), SVector(3, 1), SVector(1, 2), SVector(3, 2), SVector(4, 2)
get_edge_list(::Pyr)  = SVector(1, 2), SVector(2, 4), SVector(3, 4), SVector(3, 1), SVector(1, 5), SVector(2, 5), SVector(3, 5), SVector(4, 5)

# edges = recursive ±xyz ordering on each face
# 
#    6 ---- 8
#  / |    / |
# 5 ---- 7  |
# | 2 ___| 4
# |/     |/ 
# 1 ---- 3      
 
get_edge_list(::Hex) = SVector(1, 5), SVector(3, 7), SVector(2, 6), SVector(4, 8), 
                       SVector(1, 3), SVector(2, 4), SVector(1, 2), SVector(3, 4), 
                       SVector(5, 7), SVector(6, 8), SVector(5, 6), SVector(7, 8)


get_vertices(elem::AbstractElemShape) = equi_nodes(elem, 1)
get_vertex_fxns(elem::AbstractElemShape) =
    (rst...) -> vandermonde(elem, 1, rst...) / vandermonde(elem, 1, equi_nodes(elem, 1)...)

# assumes r1D in [-1,1], v1, v2 are vertices
function interp_1D_to_line(r1D, v1, v2)
    runit = @. (1 + r1D) / 2
    return map((a,b) -> (b-a) * runit .+ a, v1, v2)
end

"""
    interp_1D_to_edges(elem::AbstractElemShape, r1D)

Interpolates points r1D to the edges of an element (elem = :Tri, :Pyr, or :Tet)
"""
function interp_1D_to_edges(elem::AbstractElemShape, r1D)
    v = get_vertices(elem)
    edges = get_edge_list(elem)
    edge_pts = ntuple(x -> zeros(length(r1D), length(edges)), length(v))
    for (i, e) in enumerate(edges)
        edge_pts_i = interp_1D_to_line(r1D, getindex.(v, e[1]), getindex.(v, e[2]))
        setindex!.(edge_pts, edge_pts_i, :, i)
    end
    return vec.(edge_pts)
end

"""
    edge_basis(elem::AbstractElemShape, N, rst...)

Returns the generalized Vandermonde matrix evaluated using an edge basis (e.g.,
degree `N` polynomials over an edge, but linearly blended into the interior). The 
dimension of the resulting space is simply the number of total nodes on edges of 
a degree `N` element. 
"""
function edge_basis(elem::AbstractElemShape, N, rst...)
    vertices = get_vertices(elem)
    edges = get_edge_list(elem)
    vertex_functions = get_vertex_fxns(elem)
    return edge_basis(N, vertices, edges, (N, r)->basis(Line(), N, r), vertex_functions, rst...)
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

    V = zeros(length(first(rst)), length(edges) * (N-1) + length(first(vertices)))
    V[:,1:size(V1, 2)] .= V1 # initialize vertex functions

    id = size(V1,2) + 1
    for e in edges
        # edge parametrization from paper:
        # "A Comparison of High-Order Interpolation Nodes for the Pyramid"
        r1D_edge = V1[:,e[1]] - V1[:,e[2]]
        V1D, _ = basis1D(N-2, r1D_edge)
        for i in axes(V1D, 2)
            V[:, id] = V1D[:, i] .* V1[:, e[1]] .* V1[:, e[2]]
            id += 1
        end
    end
    return V
end

"""
    function find_face_nodes(elem, r, s, tol=50*eps())
    function find_face_nodes(elem, r, s, t, tol=50*eps())

Given volume nodes `r`, `s`, `t`, finds face nodes. Note that this function implicitly
defines an ordering on the faces.
"""
function find_face_nodes(::Tri, r, s, tol=50*eps())
    e1 = findall(@. abs(s + 1) < tol)
    e2 = findall(@. abs(r + s) < tol)
    e3 = findall(@. abs(r + 1) < tol)
    return e1, e2, reverse(e3)
end

function find_face_nodes(::Quad, r, s, tol=50*eps())
    e1 = findall(@. abs(r + 1) < tol)
    e2 = findall(@. abs(r - 1) < tol)
    e3 = findall(@. abs(s + 1) < tol)
    e4 = findall(@. abs(s - 1) < tol)
    return e1, e2, e3, e4
end

function find_face_nodes(::Hex, r, s, t, tol=50*eps())
    fv1 = findall(@. abs(r + 1) < tol)
    fv2 = findall(@. abs(r - 1) < tol)
    fv3 = findall(@. abs(s + 1) < tol)
    fv4 = findall(@. abs(s - 1) < tol)
    fv5 = findall(@. abs(t + 1) < tol)
    fv6 = findall(@. abs(t - 1) < tol)
    return fv1, fv2, fv3, fv4, fv5, fv6
end

function find_face_nodes(::Tet, r, s, t, tol=50*eps())
    fv1 = findall(@. abs(s + 1) < tol)
    fv2 = findall(@. abs(r + s + t + 1) < tol)
    fv3 = findall(@. abs(r + 1) < tol)
    fv4 = findall(@. abs(t + 1) < tol)
    return fv1, fv2, fv3, fv4
end

# Faces are ordered as described in "Coarse mesh partitioning for tree based AMR" 
# by Burstedde and Holke (2018). https://arxiv.org/pdf/1611.02929.pdf
function find_face_nodes(::Wedge, r, s, t, tol=50*eps())
    fv1 = findall(@. abs(s + 1) < tol)  # first quad face
    fv2 = findall(@. abs(r + s) < tol)  # second quad face
    fv3 = findall(@. abs(r + 1) < tol)  # third quad face
    fv4 = findall(@. abs(t + 1) < tol)  # bottom tri face
    fv5 = findall(@. abs(t - 1) < tol)  # top tri face
    return fv1, fv2, fv3, fv4, fv5
end

function find_face_nodes(::Pyr, r, s, t, tol=50*eps())
    fv1 = findall(@. abs(r + 1) < tol)  # +/- r tri faces
    fv2 = findall(@. abs(r + t) < tol)   
    fv3 = findall(@. abs(s + 1) < tol)  # +/- s tri faces
    fv4 = findall(@. abs(s + t) < tol)  
    fv5 = findall(@. abs(t + 1) < tol)  # bottom quad face
    return fv1, fv2, fv3, fv4, fv5
end

# face vertices = face nodes of degree 1
face_vertices(::Line) = 1, 2
face_vertices(elem) = find_face_nodes(elem, nodes(elem, 1)...)

# for 2D elements, an edge is the same as a face
face_basis(elem::T, N, rst...) where {T <: Union{Tri, Quad}} = edge_basis(elem, N, rst...)

num_vertices(::Hex) = 8
num_faces(::Hex) = 6
face_type(::Hex) = Quad()
num_face_only_nodes(::Hex, N) = (N - 1)^2

num_vertices(::Tet) = 4
num_faces(::Tet) = 4
face_type(::Tet) = Tri()
# 3 edges per face * (N-1) nodes per edge, 3 vertices
num_face_only_nodes(::Tet, N) = max(0, (N + 1) * (N + 2) ÷ 2 - (3 * (N-1) + 3)) 

function pointwise_product_of_columns(A)
    a = ones(size(A, 1))
    for A_i in eachcol(A)
        @. a *= A_i
    end
    return a
end

# 3D face basis
function face_basis(elem, N, r, s, t)

    if (elem isa Wedge) || (elem isa Pyr)
        @error "Face bases for wedges and pyramids not yet supported"
    end
        
    V_edge = edge_basis(elem, N, r, s, t)    
    if (N < 2 && elem isa Hex) || (N < 3 && elem isa Tet)
        return V_edge    
    end

    # initialize vertex and edge basis functions
    V = zeros(length(r), size(V_edge, 2) + num_faces(elem) * num_face_only_nodes(elem, N))
    V[:, 1:size(V_edge, 2)] .= V_edge
    id = size(V_edge, 2) + 1

    V1 = view(V_edge, :, 1:num_vertices(elem))
    r1, s1 = nodes(face_type(elem), 1)
    face_vertices = NodesAndModes.face_vertices(elem)
    for fids in face_vertices

        # compute face coordinates as a barycentric combo of reference vertices
        rf = sum([V1[:, fids[i]] * r1[i] for i in eachindex(r1, s1)])
        sf = sum([V1[:, fids[i]] * s1[i] for i in eachindex(r1, s1)])
        
        # extend face polynomials linearly
        linear_face_basis = pointwise_product_of_columns(V1[:, fids])
        if elem isa Hex
            V_face = vandermonde(face_type(elem), N-2, rf, sf)
        elseif elem isa Tet
            V_face = vandermonde(face_type(elem), N-3, rf, sf)
        end 
        for i in axes(V_face, 2)
            @. V[:, id] = V_face[:, i] * linear_face_basis
            id += 1
        end
    end

    return V
end

"""
    build_warped_nodes(elem::AbstractElemShape, N, r1D)

Computes degree N warp-and-blend interpolation nodes for elem = Tri(), Pyr(), or
Tet() based on the 1D node set "r1D". Returns a tuple "rst" containing arrays of
interpolation points.
"""
function build_warped_nodes(elem::AbstractElemShape, N, r1D)
    r1D_equi = equi_nodes(Line(), N)
    rst_edge_equi = interp_1D_to_edges(elem, r1D_equi)
    V_edge = edge_basis(elem, N, rst_edge_equi...)
    rst_edge = interp_1D_to_edges(elem,r1D)

    c = (x->V_edge \ x).(rst_edge) # should be lsq solve

    rst_equi = equi_nodes(elem, N)
    return (x->edge_basis(elem, N, rst_equi...) * x).(c)
end
