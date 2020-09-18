"""
    gauss_lobatto_quad(α, β, N)

Initialize Legendre-Gauss-Lobatto quadrature points.

# Examples
```jldoctest
"""
function gauss_lobatto_quad(α, β, N)

    @assert (α==0) & (β==0) "alpha/beta must be zero for Lobatto"
    # if
    #     error("alpha/beta not zero")
    # end
    x = zeros(N+1, 1)
    w = zeros(N+1, 1)
    if N == 0
        x[1] = 0
        w[1] = 2
    elseif N == 1
        x[1] = -1.0
        x[2] = 1.0
        w[1] = 1.0
        w[2] = 1.0
    else
        xint, w = gauss_quad(α+1, β+1, N-2)
        x = [-1 transpose(xint) 1]

        V = vandermonde_1D(N,x)
        w = vec(sum(inv(V*V'),dims=2))
    end
    return x[:],w[:]
end


"""
    gauss_quad(α, β, N)

Initialize weights and nodes (w,x) of Gaussian quadrature of Jacobi Polynomial
(α,β)

# Examples
```jldoctest
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
    grad_jacobiP(r, α, β, N)

Evaluate derivative of Jacobi Polynomial (α, β) of order N at r

# Examples
```jldoctest
"""

function grad_jacobiP(r, α, β, N)
    dP = zeros(length(r), 1)
    if N != 0
        dP = sqrt(N*(N+α+β+1))*jacobiP(r,α+1,β+1,N-1)
    end
    return dP
end

"""
    jacobiP(x, α, β, N)

Evaluate Jacobi Polynomial (α, β) of order N at x

# Examples
```jldoctest
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
    return P;
end


"""
    V,Vr = basis_1D(N,r)
"""
function basis_1D(N,r)
    V1D,Vr1D = ntuple(x->zeros(length(r), N+1),2)
    for j = 1:N+1
        V1D[:,j] = jacobiP(r[:], 0, 0, j-1)
        Vr1D[:,j] = grad_jacobiP(r[:], 0, 0, j-1)
    end
    return V1D,Vr1D
end

"""
    vandermonde_1D(N, r)

Initialize the 1D Vandermonde matrix of order N Legendre polynomials at nodes r

# Examples
N = 2
r,w = gauss_lobatto_quad(0,0,N)
V = vandermonde_1D(N,r)
```jldoctest
"""
function vandermonde_1D(N, r)
    V,Vr = basis_1D(N,r)
    return V
end

"""
    grad_vandermonde_1D(N, r)

Initialize the 1D Vandermonde matrix of order N Legendre polynomials at nodes r

# Examples
N = 2
r,w = gauss_lobatto_quad(0,0,N)
Vr = grad_vandermonde_1D(N,r)
```jldoctest
"""
function grad_vandermonde_1D(N, r)
    V,Vr = basis_1D(N,r)
    return Vr
end
