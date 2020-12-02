module NodesAndModes

import LinearAlgebra: diagm, eigen
import SpecialFunctions: gamma

# export some convenient 1D routines by default
include("nodes_and_modes_1D.jl")
export gauss_lobatto_quad, gauss_quad
export jacobiP, grad_jacobiP

include("warpblend_interp_nodes.jl")
include("meshgrid.jl")

# export submodules for each element type
##################
# 1D elements
##################
export Line

##################
# 2D elements
##################
export Tri
export Quad

##################
# 3D elements
##################
export Hex
export Wedge
export Pyr
export Tet

#####
##### Submodule for 1D interval
#####

module Line
import LinearAlgebra: diagm, eigen
import SpecialFunctions: gamma
using ..NodesAndModes
include("line_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for triangles
#####

module Tri
using DelimitedFiles # to read quadrature node data
using ..NodesAndModes
import ..build_warped_nodes
include("tri_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for quads
#####
module Quad
using ..NodesAndModes
import ..meshgrid
include("quad_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for tets
#####
module Tet
using DelimitedFiles # to read quadrature node data
using ..NodesAndModes
import ..meshgrid
import ..build_warped_nodes
include("tet_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for pyr
#####
module Pyr
using ..NodesAndModes
import ..meshgrid
import ..build_warped_nodes
include("pyr_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for wedge (prism)
#####
module Wedge
import ..meshgrid
using ..NodesAndModes
include("wedge_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

#####
##### Submodule for hexes
#####
module Hex
import ..meshgrid
using ..NodesAndModes
include("hex_element.jl")
include("vdm.jl") # VDM and grad_VDM
export vandermonde, grad_vandermonde, basis
export nodes, equi_nodes, quad_nodes
end

end # module
