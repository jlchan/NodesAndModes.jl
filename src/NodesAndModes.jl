module NodesAndModes

using DelimitedFiles # to read quadrature node data
import LinearAlgebra: diagm, eigen
import SpecialFunctions: gamma
include("meshgrid.jl")

# each type of element shape - used for dispatch only
abstract type AbstractElemShape end
struct Line <: AbstractElemShape end
struct Quad <: AbstractElemShape end
struct Tri <: AbstractElemShape end
struct Hex <: AbstractElemShape end
struct Wedge <: AbstractElemShape end
struct Pyr <: AbstractElemShape end
struct Tet <: AbstractElemShape end

export AbstractElemShape
export Line
export Tri, Quad
export Hex, Wedge, Pyr, Tet

# shared functions
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
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
include("tet_element.jl")

end # module
