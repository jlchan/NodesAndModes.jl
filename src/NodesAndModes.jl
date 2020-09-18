module NodesAndModes
using LinearAlgebra
using SpecialFunctions

include("nodes_and_modes_1D.jl")

# export 1D routines by default (there's only one type of element in 1D)
export gauss_lobatto_quad, gauss_quad
export jacobiP, grad_jacobiP
export vandermonde_1D, grad_vandermonde_1D
export basis_1D # returns all VDMs

#export submodules
export Tri
export Quad
export Hex
export Pyr
export Wedge

#####
##### Submodule for triangles
#####

module Tri
using DelimitedFiles # to read quadrature node data
using ..NodesAndModes
include("nodes_and_modes_2D_tri.jl")
include("warpblend_interp_nodes.jl")
export vandermonde_2D, grad_vandermonde_2D
export basis_2D
export nodes_2D, equi_nodes_2D, quad_nodes_2D
end

#####
##### Submodule for quads
#####
module Quad
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_2D_quad.jl")
export vandermonde_2D, grad_vandermonde_2D
export basis_2D
export nodes_2D, equi_nodes_2D, quad_nodes_2D
end

# #####
# ##### Submodule for tets
# #####
# module Tet
# using ..NodesAndModes
# # include("nodes_and_modes_3D_tet.jl")
# # export vandermonde_3D, grad_vandermonde_3D, basis_3D
# # export nodes_3D, equi_nodes_3D, quad_nodes_3D
# end

#####
##### Submodule for pyr
#####
module Pyr
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_pyr.jl")
export vandermonde_3D, grad_vandermonde_3D, basis_3D
export equi_nodes_3D, quad_nodes_3D
# export nodes_3D
end

#####
##### Submodule for wedge (prism)
#####
module Wedge
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_wedge.jl")
export vandermonde_3D, grad_vandermonde_3D, basis_3D
export equi_nodes_3D, quad_nodes_3D, nodes_3D
end

#####
##### Submodule for hexes
#####
module Hex
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_hex.jl")
export vandermonde_3D, grad_vandermonde_3D
export basis_3D
export nodes_3D, equi_nodes_3D, quad_nodes_3D
end

end # module
