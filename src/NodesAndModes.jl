module NodesAndModes
using LinearAlgebra
using SpecialFunctions

# export some convenient 1D routines by default
include("nodes_and_modes_1D.jl")
export gauss_lobatto_quad, gauss_quad
export jacobiP, grad_jacobiP

#export submodules
export Line # 1D
export Tri #2D
export Quad
export Hex #3D
export Wedge
export Pyr
export Tet

#####
##### Submodule for 1D interval
#####

module Line
using SpecialFunctions
using ..NodesAndModes
include("nodes_and_modes_1D.jl")
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for triangles
#####

module Tri
using DelimitedFiles # to read quadrature node data
using ..NodesAndModes
include("nodes_and_modes_2D_tri.jl")
include("warpblend_interp_nodes.jl")
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for quads
#####
module Quad
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_2D_quad.jl")
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
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
export vandermonde, grad_vandermonde, basis
export equi_nodes, quad_nodes
# export nodes
end

#####
##### Submodule for pyr
#####
module Pyr
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_pyr.jl")
include("warpblend_interp_nodes.jl")
export vandermonde, grad_vandermonde, basis
export equi_nodes, quad_nodes
# export nodes
end

#####
##### Submodule for wedge (prism)
#####
module Wedge
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_wedge.jl")
export vandermonde, grad_vandermonde, basis
export equi_nodes, quad_nodes, nodes
end

#####
##### Submodule for hexes
#####
module Hex
import VectorizedRoutines.Matlab.meshgrid
using ..NodesAndModes
include("nodes_and_modes_3D_hex.jl")
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

end # module
