module NodesAndModes

using DelimitedFiles # to read quadrature node data
import LinearAlgebra: diagm, eigen, Symmetric
import SpecialFunctions: gamma
using StaticArrays: SVector

include("meshgrid.jl")

# each type of element shape - used for dispatch only
abstract type AbstractElementShape{NDIMS} end

abstract type AbstractTensorProductElement{NDIMS} <: AbstractElementShape{NDIMS} end
struct Line <: AbstractTensorProductElement{1} end
struct Quad <: AbstractTensorProductElement{2} end
struct Hex <: AbstractTensorProductElement{3} end

abstract type AbstractSimplexElement{NDIMS} <: AbstractElementShape{NDIMS} end
struct Tri <: AbstractSimplexElement{2} end
struct Tet <: AbstractSimplexElement{3} end

# `node_ids_by_face` is an optional container for node ids of each face, since there is more than
# one face type for both Wedge and Pyramid types.
struct Wedge{T} <: AbstractElementShape{3}
    node_ids_by_face::T
end
Wedge() = Wedge(nothing)

struct Pyr{T} <: AbstractElementShape{3}
    node_ids_by_face::T
end
Pyr() = Pyr(nothing)

# TODO: deprecate AbstractElemShape{NDIMS} - it's just an alias for AbstractElementShape
const AbstractElemShape{NDIMS} = AbstractElementShape{NDIMS}
export AbstractElemShape

export AbstractElementShape, AbstractTensorProductElement, AbstractSimplexElement
export Line
export Tri, Quad
export Hex, Wedge, Pyr, Tet

# shared functions
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes, stroud_quad_nodes
include("common_functions.jl") # VDM, grad_VDM, and docs

# interp nodes for Tet(), Tri()
include("warpblend_interp_nodes.jl")

# export some convenient 1D routines by default
include("nodes_and_modes_1D.jl")
export gauss_lobatto_quad, gauss_quad
export jacobiP, grad_jacobiP

##################
# 1D elements
##################
include("line_element.jl")

# ##################
# # 2D elements
# ##################
include("tri_element.jl")
include("quad_element.jl")

##################
# 3D elements
##################
include("hex_element.jl")
include("wedge_element.jl")
include("pyr_element.jl")

export jaskowiec_sukumar_quad_nodes
include("tet_element.jl")

end # module
