"""
    gauss_lobatto_quad(α, β, N)

Computes Legendre-Gauss-Lobatto quadrature points and weights with Jacobi weights α,β.
"""
function gauss_lobatto_quad(α::Int, β::Int, N)

    if !(iszero(α) && iszero(β))
        error("alpha/beta must be zero for Lobatto")
    end

    x = zeros(N+1)
    w = zeros(N+1)
    if N == 0
        x[1] = 0.
        w[1] = 2.
    elseif N == 1
        x[1] = -1.0
        x[2] = 1.0
        w[1] = 1.0
        w[2] = 1.0
    elseif N==2
        x[1] = -1.
        x[2] = 0.
        x[3] = 1.
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

Compute nodes and weights for Gaussian quadrature with Jacobi weights (α,β).
"""
function gauss_quad(α, β, N)
    if N == 0
        x = [-(α-β)/(α+β+2)]
        w = [2.]
        return x, w
    end

    J = zeros(N+1, N+1)
    h₁ = @. 2*(0:N)+α+β    
    J = diagm(0 => @. -1/2*(α^2-β^2)/(h₁+2)/h₁) + 
        diagm(1 => @. 2/(h₁[1:N]+2) * sqrt((1:N) * ((1:N)+α+β) * ((1:N)+α) * ((1:N)+β) / (h₁[1:N]+1) / (h₁[1:N]+3)))
    if α + β < 10*eps()
        J[1,1] = 0.0
    end

    x, V = eigen(Symmetric(J+J'))
    w = @. transpose(V[1,:])^2 * 2^(α+β+1) / (α+β+1) * gamma(α+1) * gamma(β+1) / gamma(α+β+1)
    return x[:], w[:]
end

"""
    jacobiP(x, α, β, N)

Evaluate the Jacobi Polynomial (α, β) of order N at points x
"""
function jacobiP(x, α, β, N)

    PL = zeros(N+1, length(x))
    γ₀ = 2^(α+β+1) / (α+β+1) * gamma(α+1) * gamma(β+1) / gamma(α+β+1)
    @. PL[1,:] = 1.0/sqrt(γ₀)
    if N == 0
        return vec(PL)
    end

    γ₁ = (α+1) * (β+1) / (α+β+3) * γ₀
    @. PL[2, :] = ((α+β+2) * x/2 + (α-β)/2) / sqrt(γ₁)
    if N == 1        
        return PL[N+1,:]
    end

    aold = 2 / (2+α+β) * sqrt((α+1) * (β+1) / (α + β + 3))

    for i = 1:N-1
        h₁ = 2i+α+β
        anew = 2 / (h₁+2) * sqrt((i+1) * (i+1+α+β) * (i+1+α) * (i+1+β) / (h₁+1) / (h₁+3))
        bnew = -(α^2-β^2) / h₁ / (h₁+2)
        @. PL[i+2,:] = 1 / anew * (-aold * PL[i,:] .+ (x .- bnew) .* PL[i+1,:])
        aold = anew
    end
    
    return PL[N+1,:]
end

"""
    grad_jacobiP(r, α, β, N)

Evaluate derivative of Jacobi Polynomial (α, β) of order N at points r.
"""

function grad_jacobiP(r, α, β, N)    
    if N != 0
        return sqrt(N * (N+α+β+1)) * jacobiP(r, α+1, β+1, N-1)
    else
        return zeros(length(r))
    end
end
