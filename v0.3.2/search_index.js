var documenterSearchIndex = {"docs":
[{"location":"Tri/#NodesAndModes.Tri","page":"Triangular element","title":"NodesAndModes.Tri","text":"","category":"section"},{"location":"Tri/","page":"Triangular element","title":"Triangular element","text":"CurrentModule = NodesAndModes.Tri","category":"page"},{"location":"Tri/","page":"Triangular element","title":"Triangular element","text":"Modules = [Tri]","category":"page"},{"location":"Tri/#NodesAndModes.Tri.equi_nodes-Tuple{Any}","page":"Triangular element","title":"NodesAndModes.Tri.equi_nodes","text":"equi_nodes(N)\n\nComputes equally spaced nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.grad_vandermonde-Tuple{Any,Any}","page":"Triangular element","title":"NodesAndModes.Tri.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Triangular element","title":"NodesAndModes.Tri.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.nodes-Tuple{Any}","page":"Triangular element","title":"NodesAndModes.Tri.nodes","text":"nodes(N)\n\nComputes interpolation nodes of degree N. Edge nodes coincide with (N+1)-point Lobatto points.\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.rstoab-Tuple{Any,Any}","page":"Triangular element","title":"NodesAndModes.Tri.rstoab","text":"rstoab(r, s)\n\nConverts from reference bi-unit right triangle coordinate (r,s) to polynomial basis evaluation coordinates (a,b) on the domain [-1,1]^2\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.simplex_2D-NTuple{4,Any}","page":"Triangular element","title":"NodesAndModes.Tri.simplex_2D","text":"simplex_2D(a, b, i, j)\n\nEvaluate 2D PKDO basis phi_ij at points (a,b) on the Duffy domain [-1,1]^2\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.stroud_quad_nodes-Tuple{Any}","page":"Triangular element","title":"NodesAndModes.Tri.stroud_quad_nodes","text":"stroud_quad_nodes(N)\n\nReturns Stroud-type quadrature nodes constructed from the tensor product of (N+1)-point Gauss-Jacobi rules. Exact for degree 2N polynomials\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.vandermonde-Tuple{Any,Any}","page":"Triangular element","title":"NodesAndModes.Tri.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Tri/#NodesAndModes.Tri.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Triangular element","title":"NodesAndModes.Tri.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes","page":"Helper functions","title":"NodesAndModes","text":"","category":"section"},{"location":"nodes_and_modes_1D/","page":"Helper functions","title":"Helper functions","text":"CurrentModule = NodesAndModes","category":"page"},{"location":"nodes_and_modes_1D/","page":"Helper functions","title":"Helper functions","text":"Modules = [NodesAndModes]","category":"page"},{"location":"nodes_and_modes_1D/#NodesAndModes.build_warped_nodes-Tuple{Any,Symbol,Any}","page":"Helper functions","title":"NodesAndModes.build_warped_nodes","text":"build_warped_nodes(N,elem::Symbol,r1D)\n\nComputes degree N warp-and-blend interpolation nodes for elem = :Tri, :Pyr, or :Tet based on the 1D node set \"r1D\". Returns a tuple \"rst\" containing arrays of interpolation points.\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.edge_basis-Tuple{Any,Any,Any,Any,Any,Vararg{Any,N} where N}","page":"Helper functions","title":"NodesAndModes.edge_basis","text":"edge_basis(N, vertices, edges, basis, vertex_functions, rst...)\n\nComputes edge basis given vertex functions and 1D basis.\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.edge_basis-Tuple{Any,Symbol,Vararg{Any,N} where N}","page":"Helper functions","title":"NodesAndModes.edge_basis","text":"edge_basis(N, elem::Symbol, rst...)\n\nreturns generalized Vandermonde matrix evaluated using an edge basis.\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.gauss_lobatto_quad-Tuple{Any,Any,Any}","page":"Helper functions","title":"NodesAndModes.gauss_lobatto_quad","text":"gauss_lobatto_quad(α, β, N)\n\nComputes Legendre-Gauss-Lobatto quadrature points and weights with Jacobi weights α,β.\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.gauss_quad-Tuple{Any,Any,Any}","page":"Helper functions","title":"NodesAndModes.gauss_quad","text":"gauss_quad(α, β, N)\n\nCompute nodes and weights for Gaussian quadrature with Jacobi weights (α,β).\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.get_edge_list-Tuple{Symbol}","page":"Helper functions","title":"NodesAndModes.get_edge_list","text":"get_edge_list(elem::Symbol)\n\nReturns list of edges for a specific element (elem = :Tri, :Pyr, or :Tet).\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.interp_1D_to_edges-Tuple{Any,Symbol}","page":"Helper functions","title":"NodesAndModes.interp_1D_to_edges","text":"interp_1D_to_edges(r1D,elem::Symbol)\n\nInterpolates points r1D to the edges of an element (elem = :Tri, :Pyr, or :Tet)\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.jacobiP-NTuple{4,Any}","page":"Helper functions","title":"NodesAndModes.jacobiP","text":"jacobiP(x, α, β, N)\n\nEvaluate the Jacobi Polynomial (α, β) of order N at points x\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.meshgrid-Tuple{AbstractArray{T,1} where T}","page":"Helper functions","title":"NodesAndModes.meshgrid","text":"meshgrid(vx) Computes an (x,y)-grid from the vectors (vx,vx). For more information, see the MATLAB documentation.\n\nCopied and pasted directly from VectorizedRoutines.jl. Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.meshgrid-Union{Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1},AbstractArray{T,1}}} where T","page":"Helper functions","title":"NodesAndModes.meshgrid","text":"meshgrid(vx,vy,vz) Computes an (x,y,z)-grid from the vectors (vx,vy,vz). For more information, see the MATLAB documentation.\n\nCopied and pasted directly from VectorizedRoutines.jl. Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl\n\n\n\n\n\n","category":"method"},{"location":"nodes_and_modes_1D/#NodesAndModes.meshgrid-Union{Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1}}} where T","page":"Helper functions","title":"NodesAndModes.meshgrid","text":"meshgrid(vx,vy) Computes an (x,y)-grid from the vectors (vx,vy). For more information, see the MATLAB documentation.\n\nCopied and pasted directly from VectorizedRoutines.jl. Using VectorizedRoutines.jl directly causes Pkg versioning issues with SpecialFunctions.jl\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex","page":"Hexahedra element","title":"NodesAndModes.Hex","text":"","category":"section"},{"location":"Hex/","page":"Hexahedra element","title":"Hexahedra element","text":"CurrentModule = NodesAndModes.Hex","category":"page"},{"location":"Hex/","page":"Hexahedra element","title":"Hexahedra element","text":"Modules = [Hex]","category":"page"},{"location":"Hex/#NodesAndModes.Hex.basis-NTuple{4,Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.basis","text":"basis(N, r, s, t)\n\nComputes orthonormal basis of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.equi_nodes-Tuple{Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.equi_nodes","text":"equi_nodes(N)\n\nCompute equispaced nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.grad_vandermonde-Tuple{Any,Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Hexahedra element","title":"NodesAndModes.Hex.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.nodes-Tuple{Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.nodes","text":"nodes(N)\n\nComputes optimized interpolation nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.quad_nodes-Tuple{Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.quad_nodes","text":"quad_nodes(N)\n\nCompute quadrature nodes and weights exact for degree 2N+1 polynomials. Uses a tensor product Gauss quadrature rule.\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.vandermonde-Tuple{Any,Any}","page":"Hexahedra element","title":"NodesAndModes.Hex.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Hex/#NodesAndModes.Hex.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Hexahedra element","title":"NodesAndModes.Hex.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr","page":"Pyramid element","title":"NodesAndModes.Pyr","text":"","category":"section"},{"location":"Pyr/","page":"Pyramid element","title":"Pyramid element","text":"CurrentModule = NodesAndModes.Pyr","category":"page"},{"location":"Pyr/","page":"Pyramid element","title":"Pyramid element","text":"Modules = [Pyr]","category":"page"},{"location":"Pyr/#NodesAndModes.Pyr.abctorst-Tuple{Any,Any,Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.abctorst","text":"abctorst(a,b,c)\n\nConverts from Stroud coordinates (a,b,c) on [-1,1]^3 to reference element coordinates (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.basis","page":"Pyramid element","title":"NodesAndModes.Pyr.basis","text":"basis(N,r,s,t,tol=1e-12)\n\nComputes orthonormal semi-nodal basis on the biunit pyramid element.\n\nWarning: nodal derivative matrices may contain errors for nodes at t = 1. A way to avoid this is to use weak differentiation matrices computed using quadrature rules with only interior nodes.\n\n\n\n\n\n","category":"function"},{"location":"Pyr/#NodesAndModes.Pyr.equi_nodes-Tuple{Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.equi_nodes","text":"equi_nodes(N)\n\nComputes equispaced nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.grad_vandermonde-Tuple{Any,Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Pyramid element","title":"NodesAndModes.Pyr.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.nodes-Tuple{Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.nodes","text":"nodes(N)\n\nComputes interpolation nodes of degree N. Edge nodes coincide with (N+1)-point Lobatto points. Triangular face nodes coincide with Tri.nodes(N), quadrilateral face nodes coincide with tensor product (N+1)-point Lobatto points.\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.quad_nodes-Tuple{Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.quad_nodes","text":"quad_nodes(N)\n\nComputes quadrature nodes and weights which are exact for degree 2N polynomials.\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.stroud_quad_nodes-Tuple{Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.stroud_quad_nodes","text":"stroud_quad_nodes(N)\n\nReturns Stroud-type quadrature nodes constructed from the tensor product of (N+1)-point Gauss-Jacobi rules. Exact for degree 2N polynomials\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.vandermonde-Tuple{Any,Any}","page":"Pyramid element","title":"NodesAndModes.Pyr.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Pyr/#NodesAndModes.Pyr.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Pyramid element","title":"NodesAndModes.Pyr.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge","page":"Wedge element","title":"NodesAndModes.Wedge","text":"","category":"section"},{"location":"Wedge/","page":"Wedge element","title":"Wedge element","text":"CurrentModule = NodesAndModes.Wedge","category":"page"},{"location":"Wedge/","page":"Wedge element","title":"Wedge element","text":"Modules = [Wedge]","category":"page"},{"location":"Wedge/#NodesAndModes.Wedge.basis","page":"Wedge element","title":"NodesAndModes.Wedge.basis","text":"basis(N, r, s, t)\n\nComputes orthonormal basis of degree N at points (r,s,t)\n\n\n\n\n\n","category":"function"},{"location":"Wedge/#NodesAndModes.Wedge.equi_nodes-Tuple{Any}","page":"Wedge element","title":"NodesAndModes.Wedge.equi_nodes","text":"equi_nodes(N)\n\nCompute equispaced nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.grad_vandermonde-Tuple{Any,Any}","page":"Wedge element","title":"NodesAndModes.Wedge.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Wedge element","title":"NodesAndModes.Wedge.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.nodes-Tuple{Any}","page":"Wedge element","title":"NodesAndModes.Wedge.nodes","text":"nodes(N)\n\nComputes interpolation nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.quad_nodes-Tuple{Any}","page":"Wedge element","title":"NodesAndModes.Wedge.quad_nodes","text":"quad_nodes(N)\n\nReturns quadrature nodes and weights which exactly integrate degree 2N polynomials.\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.vandermonde-Tuple{Any,Any}","page":"Wedge element","title":"NodesAndModes.Wedge.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Wedge/#NodesAndModes.Wedge.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Wedge element","title":"NodesAndModes.Wedge.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"authors/#Authors","page":"Authors","title":"Authors","text":"","category":"section"},{"location":"authors/","page":"Authors","title":"Authors","text":"The original port from Matlab to Julia was done by Yimin Lin. Subsequent development was done by Jesse Chan.","category":"page"},{"location":"Tet/#NodesAndModes.Tet","page":"Tetrahedral element","title":"NodesAndModes.Tet","text":"","category":"section"},{"location":"Tet/","page":"Tetrahedral element","title":"Tetrahedral element","text":"CurrentModule = NodesAndModes.Tet","category":"page"},{"location":"Tet/","page":"Tetrahedral element","title":"Tetrahedral element","text":"Modules = [Tet]","category":"page"},{"location":"Tet/#NodesAndModes.Tet.basis-NTuple{4,Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.basis","text":"basis(N, r, s, t)\n\nComputes orthonormal basis of degree N at points (r,s,t)\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.grad_vandermonde-Tuple{Any,Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Tetrahedral element","title":"NodesAndModes.Tet.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.nodes-Tuple{Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.nodes","text":"nodes(N)\n\nComputes interpolation nodes of degree N.  Edge nodes coincide with (N+1)-point Lobatto points, while face nodes coincide with Tr.nodes(N).\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.quad_nodes_tet-Tuple{Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.quad_nodes_tet","text":"quad_nodes_tet(N)\n\nReturns quadrature nodes and weights which exactly integrate degree N polynomials\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.simplex_3D-NTuple{6,Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.simplex_3D","text":"simplex_3D(a, b, c, i, j, k)\n\nEvaluate 3D \"Legendre\" basis phi_ijk at (a,b,c) coordinates on the [-1,1] cube\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.stroud_quad_nodes-Tuple{Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.stroud_quad_nodes","text":"stroud_quad_nodes(N)\n\nReturns Stroud-type quadrature nodes constructed from the tensor product of (N+1)-point Gauss-Jacobi rules. Exact for degree 2N polynomials\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.vandermonde-Tuple{Any,Any}","page":"Tetrahedral element","title":"NodesAndModes.Tet.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Tet/#NodesAndModes.Tet.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Tetrahedral element","title":"NodesAndModes.Tet.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad","page":"Quadrilateral element","title":"NodesAndModes.Quad","text":"","category":"section"},{"location":"Quad/","page":"Quadrilateral element","title":"Quadrilateral element","text":"CurrentModule = NodesAndModes.Quad","category":"page"},{"location":"Quad/","page":"Quadrilateral element","title":"Quadrilateral element","text":"Modules = [Quad]","category":"page"},{"location":"Quad/#NodesAndModes.Quad.basis-Tuple{Any,Any,Any}","page":"Quadrilateral element","title":"NodesAndModes.Quad.basis","text":"basis(N,r,s)\n\nComputes orthonormal basis of degree N at coordinates (r,s)\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad.grad_vandermonde-Tuple{Any,Any}","page":"Quadrilateral element","title":"NodesAndModes.Quad.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Quadrilateral element","title":"NodesAndModes.Quad.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad.quad_nodes-Tuple{Any}","page":"Quadrilateral element","title":"NodesAndModes.Quad.quad_nodes","text":"quad_nodes(N)\n\nCompute quadrature nodes and weights exact for degree 2N+1 polynomials. Uses a tensor product Gauss quadrature rule.\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad.vandermonde-Tuple{Any,Any}","page":"Quadrilateral element","title":"NodesAndModes.Quad.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Quad/#NodesAndModes.Quad.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Quadrilateral element","title":"NodesAndModes.Quad.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"#NodesAndModes","page":"Home","title":"NodesAndModes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NodesAndModes.jl is a package to compute nodes (interpolation and quadrature points) and modes (orthogonal polynomials) on various reference elements for use in high order finite element and nodal discontinuous Galerkin (DG) methods.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The codes inspired by the Matlab codes for the book \"Nodal Discontinuous Galerkin methods\" by Hesthaven and Warburton (2007) and high order interpolation nodes on triangles, tetrahedra, and pyramids using the \"Interpolatory Warp and Blend\" procedure from Chan and Warburton 2015.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is registered and can be installed via ] add NodesAndModes. Julia v1.4 is required.","category":"page"},{"location":"#Package-organization","page":"Home","title":"Package organization","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NodesAndModes.jl supports seven types of elements in 1D, 2D, and 3D. Routines for each element are packaged into the following submodules","category":"page"},{"location":"","page":"Home","title":"Home","text":"NodesAndModes.Line (1D lines/intervals)\nNodesAndModes.Tri (2D triangles)\nNodesAndModes.Quad (2D quadrilaterals)\nNodesAndModes.Tet (3D tetrahedra)\nNodesAndModes.Hex (3D hexahedra)\nNodesAndModes.Wedge (3D wedges/prisms)\nNodesAndModes.Pyr (3D pyramids)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Each module exports the following functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"basis: computes Vandermonde matrix (columns are evaluations of orthogonal polynomials at different points) and derivative matrices (columns are derivatives of orthogonal polynomials at different points) constructed using orthogonal polynomials on the reference element\nnodes: computes (non-uniform) interpolation nodes on the reference element. All interpolation nodes\nquad_nodes: computes quadrature nodes and weights on the reference element. quad_nodes(N) returns a quadrature rule exact for degree 2N polynomials (e.g., exact integration of the mass matrix).\nequi_nodes: computes equispaced nodes on the reference element (for plotting)","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nEach module also exports vandermonde and grad_vandermonde for similarity with the Matlab codes in Hesthaven/Warburton 2007. These just call the basis routine and return either Vandermonde or derivative matrices.","category":"page"},{"location":"#Example-usage","page":"Home","title":"Example usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To compute a 1D Vandermonde matrix using Gauss-Lobatto points and orthonormal polynomials.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using NodesAndModes.Line\njulia> N = 2\njulia> r = nodes(N)\njulia> V = vandermonde(N,r)","category":"page"},{"location":"","page":"Home","title":"Home","text":"To compute a 2D triangular Vandermonde matrix from Warp & Blend points (see Warburton 2006) and orthonormal polynomials on the triangle (with coordinates r,s)","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using NodesAndModes.Tri\njulia> N = 2\njulia> r,s = nodes(N)\njulia> V = vandermonde(N,r,s)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Nodal differentiation matrices Dr and Ds can be computed via","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using NodesAndModes.Tri\njulia> N = 2\njulia> r,s = nodes(N)\njulia> V,Vr,Vs = basis(N,r,s)\njulia> Dr,Ds = (A->A/V).((Vr,Vs))","category":"page"},{"location":"","page":"Home","title":"Home","text":"such that Dr*f(r,s) ≈ df/dr. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"A mass matrix M and weak differentation matrices Qr,Qs in finite element or DG methods can be computed using quadrature via","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using LinearAlgebra\njulia> using NodesAndModes.Tri\njulia> N = 2\njulia> r,s = nodes(N)\njulia> V = vandermonde(N,r,s)\njulia> rq,sq,wq = quad_nodes(N)\njulia> Vq,Vrq,Vsq = (A->A/V).(basis(N,rq,sq))\njulia> M,Qr,Qs = (A->Vq'*diagm(wq)*A).((Vq,Vrq,Vsq))","category":"page"},{"location":"Line/#NodesAndModes.Line","page":"Line (1D) element","title":"NodesAndModes.Line","text":"","category":"section"},{"location":"Line/","page":"Line (1D) element","title":"Line (1D) element","text":"CurrentModule = NodesAndModes.Line","category":"page"},{"location":"Line/","page":"Line (1D) element","title":"Line (1D) element","text":"Modules = [Line]","category":"page"},{"location":"Line/#NodesAndModes.Line.basis-Tuple{Any,Any}","page":"Line (1D) element","title":"NodesAndModes.Line.basis","text":"basis(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N (along with the derivative matrix Vr) at points r.\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.equi_nodes-Tuple{Any}","page":"Line (1D) element","title":"NodesAndModes.Line.equi_nodes","text":"equi_nodes(N)\n\nComputes equally spaced nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.grad_vandermonde-Tuple{Any,Any}","page":"Line (1D) element","title":"NodesAndModes.Line.grad_vandermonde","text":"grad_vandermonde(N,r)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.grad_vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Line (1D) element","title":"NodesAndModes.Line.grad_vandermonde","text":"grad_vandermonde(N,rst...)\n\nComputes the generalized Vandermonde derivative matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.nodes-Tuple{Any}","page":"Line (1D) element","title":"NodesAndModes.Line.nodes","text":"nodes(N)\n\nComputes interpolation nodes of degree N.\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.quad_nodes-Tuple{Any}","page":"Line (1D) element","title":"NodesAndModes.Line.quad_nodes","text":"quad_nodes(N)\n\nComputes (N+1)-point Gauss quadrature rule (exact for degree 2N+1 polynomials)\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.vandermonde-Tuple{Any,Any}","page":"Line (1D) element","title":"NodesAndModes.Line.vandermonde","text":"vandermonde(N,r)\n\nComputes the generalized Vandermonde matrix V of degree N at point r. Specialized for the 1D case\n\n\n\n\n\n","category":"method"},{"location":"Line/#NodesAndModes.Line.vandermonde-Tuple{Any,Vararg{Any,N} where N}","page":"Line (1D) element","title":"NodesAndModes.Line.vandermonde","text":"vandermonde(N,rst...)\n\nComputes the generalized Vandermonde matrix V of degree N at points (r,s,t).\n\n\n\n\n\n","category":"method"}]
}
