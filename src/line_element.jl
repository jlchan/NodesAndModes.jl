"""
    basis(elem::Line,N,r)

Computes the generalized Vandermonde matrix V of degree N (along with the
derivative matrix Vr) at points r.
"""
function basis(elem::Line, N, r)
    V1D, Vr1D = ntuple(x->zeros(length(r), N+1),2)
    for j = 1:N+1
        V1D[:,j] .= jacobiP(r[:], 0, 0, j-1)
        Vr1D[:,j] .= grad_jacobiP(r[:], 0, 0, j-1)
    end
    return V1D, Vr1D
end

"""
    nodes(elem::Line,N)

Computes interpolation nodes of degree N.
"""
nodes(elem::Line, N) = first(gauss_lobatto_quad(0, 0, N))

"""
    equi_nodes(elem::Line,N)

Computes equally spaced nodes of degree N.
"""
equi_nodes(elem::Line, N) = collect(LinRange(-1, 1, N+1))

"""
    quad_nodes(elem::Line,N)

Computes (N+1)-point Gauss quadrature rule (exact for degree 2N+1 polynomials)

"""
quad_nodes(elem::Line, N) = gauss_quad(0, 0, N)
