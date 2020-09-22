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
# export Line # 1D - TODO: add
export Tri #2D
export Quad
export Hex #3D
export Wedge
export Pyr
export Tet

# #####
# ##### Submodule for 1D interval
# #####
#
# module Line
# using ..NodesAndModes
# export vandermonde, grad_vandermonde, basis
# export nodes, equi_nodes, quad_nodes
# end

#####
##### Submodule for triangles
#####

module Tri
using DelimitedFiles # to read quadrature node data
using ..NodesAndModes
include("nodes_and_modes_2D_tri.jl")
include("warpblend_interp_nodes.jl")
export vandermonde_2D, grad_vandermonde_2D, basis_2D
export nodes_2D, equi_nodes_2D, quad_nodes_2D
end

#####
##### Submodule for quads
#####
module Quad
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_2D_quad.jl")
export vandermonde_2D, grad_vandermonde_2D, basis_2D
export nodes_2D, equi_nodes_2D, quad_nodes_2D
end

#####
##### Submodule for tets
#####
module Tet
using DelimitedFiles # to read quadrature node data
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_tet.jl")
include("warpblend_interp_nodes.jl")
export vandermonde_3D, grad_vandermonde_3D, basis_3D
export equi_nodes_3D, quad_nodes_3D
# export nodes_3D
end

#####
##### Submodule for pyr
#####
module Pyr
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_pyr.jl")
include("warpblend_interp_nodes.jl")
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
export vandermonde_3D, grad_vandermonde_3D, basis_3D
export nodes_3D, equi_nodes_3D, quad_nodes_3D
end

end # module
