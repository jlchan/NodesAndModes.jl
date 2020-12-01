"""
    gauss_lobatto_quad(α, β, N)

Computes Legendre-Gauss-Lobatto quadrature points and weights.

# Examples
```jldoctest
julia> r,w = gauss_lobatto_quad(0,0,1)
 ([-1.0, 1.0], [1.0, 1.0])
```
"""
function gauss_lobatto_quad(α, β, N)

    @assert (α==0) & (β==0) "alpha/beta must be zero for Lobatto"
    x = zeros(N+1)
    w = zeros(N+1)
    if N == 0
        x[1] = 0
        w[1] = 2
    elseif N == 1
        x[1] = -1.0
        x[2] = 1.0
        w[1] = 1.0
        w[2] = 1.0
    elseif N==2
        x[1] = -1
        x[2] = 0
        x[3] = 1
        w[1] = 0.333333333333333
        w[2] = 1.333333333333333
        w[3] = 0.333333333333333
    else
        xint, wint = gauss_quad(α+1, β+1, N-2)
        x = vcat(-1,xint,1)
        wend = 2/(N*(N+1))
        w = vcat(wend, (@. wint / (1-xint^2)) ,wend)
    end
    return x[:],w[:]
end


"""
    gauss_quad(α, β, N)

Compute nodes and weights for Gaussian quadrature with Jacobi weights (α,β)

# Examples
```jldoctest
julia> r,w = gauss_quad(0,0,1)
 ([-0.5773502691896257, 0.5773502691896257], [0.9999999999999998, 0.9999999999999998])
```
"""
function gauss_quad(α, β, N)
    if N == 0
        x = [-(α-β)/(α+β+2)]
        w = [2]
        return x, w
    end

    J = zeros(N+1, N+1)
    h₁ = @. 2*(0:N)+α+β
    J = diagm(0 => @. -1/2*(α^2-β^2)/(h₁+2)/h₁) + diagm(1 => @. 2/(h₁[1:N]+2)*sqrt((1:N)*((1:N)+α+β)*((1:N)+α)*((1:N)+β)/(h₁[1:N]+1)/(h₁[1:N]+3)))
    if α+β<10*eps()
        J[1,1] = 0.0
    end
    J = J + transpose(J)

    x, V = eigen(J)
    w = @. transpose(V[1,:])^2*2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+1)
    return x[:], w[:]
end

"""
    jacobiP(x, α, β, N)

Evaluate the Jacobi Polynomial (α, β) of order N at points x

# Examples
```jldoctest
julia> r = [-1, 1]; jacobiP(r,0,0,1)
 2-element Array{Float64,1}:
 -1.224744871391589
  1.224744871391589
```
"""
function jacobiP(x, α, β, N)
    xp = x
    if size(xp, 2) == 1
        xp = transpose(xp)
    end

    PL = zeros(N+1,length(xp))
    γ₀ = 2^(α+β+1)/(α+β+1)*gamma(α+1)*gamma(β+1)/gamma(α+β+1)
    PL[1,:] .= 1.0/sqrt(γ₀)
    if N == 0
        P = transpose(PL)
        return P
    end
    γ₁ = (α+1)*(β+1)/(α+β+3)*γ₀
    PL[2,:] = ((α+β+2).*xp/2 .+ (α-β)/2)/sqrt(γ₁)
    if N == 1
        P = PL[N+1,:]
        return P
    end

    aold = 2/(2+α+β)*sqrt((α+1)*(β+1)/(α+β+3))

    for i = 1:N-1
        h₁ = 2i+α+β
        anew = 2/(h₁+2)*sqrt((i+1)*(i+1+α+β)*(i+1+α)*(i+1+β)/(h₁+1)/(h₁+3))
        bnew = -(α^2-β^2)/h₁/(h₁+2)
        # TODO: transpose?
        PL[i+2,:] = 1/anew*(-aold*transpose(PL[i,:]).+(xp.-bnew).*transpose(PL[i+1,:]))
        aold = anew
    end

    P = PL[N+1,:]
    return P
end

"""
    grad_jacobiP(r, α, β, N)

Evaluate derivative of Jacobi Polynomial (α, β) of order N at r

# Examples
```jldoctest
julia> r = [-1, 1]; grad_jacobiP(r,0,0,1)
 2×1 Array{Float64,2}:
 1.2247448713915892
 1.2247448713915892
```
"""

function grad_jacobiP(r, α, β, N)
    dP = zeros(length(r))
    if N != 0
        dP = sqrt(N*(N+α+β+1))*jacobiP(r,α+1,β+1,N-1)
    end
    return dP
end

"""
    basis(N,r)

Computes the generalized Vandermonde matrix V of degree N (along with the
derivative matrix Vr) at points r.

# Examples
```jldoctest
julia> N = 1; r = nodes(N);

julia> V,Vr = basis(N,r);

julia> V
 2×2 Array{Float64,2}:
 0.707107  -1.22474
 0.707107   1.22474

julia> Vr
 2×2 Array{Float64,2}:
 0.0  1.22474
 0.0  1.22474
```
"""
function basis(N,r)
    V1D,Vr1D = ntuple(x->zeros(length(r), N+1),2)
    for j = 1:N+1
        V1D[:,j] = jacobiP(r[:], 0, 0, j-1)
        Vr1D[:,j] = grad_jacobiP(r[:], 0, 0, j-1)
    end
    return V1D,Vr1D
end

"""
    nodes(N)

Computes interpolation nodes of degree N.

# Examples
```jldoctest
julia> r = nodes(3)
 4-element Array{Float64,1}:
 -1.0
 -0.4472135954999579
  0.4472135954999579
  1.0
```
"""
nodes(N) = first(gauss_lobatto_quad(0,0,N))

"""
    equi_nodes(N)

Computes equally spaced nodes of degree N.

# Examples
```jldoctest
julia> equi_nodes(2)
 3-element Array{Float64,1}:
 -1.0
  0.0
  1.0
```
"""
equi_nodes(N) = collect(LinRange(-1,1,N+1))

"""
    quad_nodes(N)

Computes (N+1)-point Gauss quadrature rule (exact for degree 2N+1 polynomials)

# Examples
```jldoctest
julia> r,w = quad_nodes(2);

julia> r
 3-element Array{Float64,1}:
 -0.7745966692414822
  7.771561172376096e-16
  0.7745966692414833

julia> w
 3-element Array{Float64,1}:
 0.555555555555557
 0.8888888888888867
 0.5555555555555558
```
"""
quad_nodes(N) = gauss_quad(0,0,N)
