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
        x = [-(α - β) / (α + β + 2)]
        w = [2.]
        return x, w
    end

    J = zeros(N+1, N+1)
    h₁ = @. 2 * (0:N) + α + β    
    J = diagm(0 => @. -1/2 * (α^2 - β^2) / (h₁ + 2) / h₁) + 
        diagm(1 => @. 2 / (h₁[1:N] + 2) * sqrt((1:N) * ((1:N) + α + β) * ((1:N) + α) * ((1:N) + β) / (h₁[1:N] + 1) / (h₁[1:N] + 3)))
    if α + β < 10 * eps()
        J[1,1] = 0.0
    end

    x, V = eigen(Symmetric(J+J'))
    w = @. transpose(V[1,:])^2 * 2^(α + β + 1) / (α + β + 1) * gamma(α + 1) * gamma(β + 1) / gamma(α + β + 1)
    return x[:], w[:]
end

"""
    jacobiP(x, α, β, N)
    jacobiP!(out, x, α, β, N)    
    jacobiP(x::T, α, β, N) where {T <: Real}

Evaluate the Jacobi Polynomial (α, β) of order N at points x. 
Can optionally specify the output array `out` to avoid allocations.
"""
# vectorized version. 
function jacobiP(x, α, β, N)
    PL = zeros(length(x))
    return jacobiP!(PL, x, α, β, N)
end

# non-allocating vector version
function jacobiP!(out, x, α, β, N)    
    for i in eachindex(x)
        out[i] = jacobiP(x[i], α, β, N)
    end
    return out
end

# non-allocating scalar version of jacobiP
function jacobiP(x::T, α, β, N) where {T <: Real}

    γ₀ = 2^(α + β + 1) / (α + β + 1) * gamma(α + 1) * gamma(β + 1) / gamma(α + β + 1)
    PL = 1.0 / sqrt(γ₀)
    if N == 0
        return PL
    end

    γ₁ = (α + 1) * (β + 1) / (α + β + 3) * γ₀
    PL_prev = PL
    PL = ((α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ₁)
    if N == 1        
        return PL
    end

    aold = 2 / (2+α+β) * sqrt((α+1) * (β+1) / (α + β + 3))
    PL_prev2 = PL_prev
    PL_prev = PL
    for i = 1:N-1
        h₁ = 2i + α + β
        anew = 2 / (h₁+2) * sqrt((i+1) * (i+1+α+β) * (i+1+α) * (i+1+β) / (h₁+1) / (h₁+3))
        bnew = -(α^2-β^2) / h₁ / (h₁+2)
        PL = 1 / anew * (-aold * PL_prev2 + (x - bnew) * PL_prev)
        aold = anew

        # reset recursion
        PL_prev2 = PL_prev
        PL_prev = PL
    end
    
    return PL
end

"""
    grad_jacobiP(r, α, β, N, tmp_array=())

Evaluate derivative of Jacobi polynomial (α, β) of order N at points r.
"""
function grad_jacobiP(r, α, β, N)
    out = zeros(length(r))
    return grad_jacobiP!(out, r, α, β, N)
end

function grad_jacobiP!(out, r, α, β, N)
    for i in eachindex(r)
        out[i] = grad_jacobiP(r[i], α, β, N)
    end
    return out
end

function grad_jacobiP(r::T, α, β, N) where {T <: Real}
    if N != 0
        return sqrt(N * (N+α+β+1)) * jacobiP(r, α+1, β+1, N-1)
    else
        return zero(typeof(r))
    end
end
