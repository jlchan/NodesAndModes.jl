"""
    gauss_lobatto_quad(α, β, N)

Initialize Legendre-Gauss-Lobatto quadrature points.

# Examples
```jldoctest
"""
function gauss_lobatto_quad(α, β, N)

    if (α!=0) & (β!=0)
        error("alpha/beta not zero")
    end
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


###
### 2D modes on triangles
###

"""
    simplex_2D(a, b, i, j)

Evaluate 2D "Legendre" basis phi_ij at (a,b) returned by rstoab

# Examples
```jldoctest
"""
function simplex_2D(a, b, i, j)
    h₁ = jacobiP(a, 0, 0, i)
    h₂ = jacobiP(b, 2i+1, 0, j)
    return @. sqrt(2.0)*h₁*h₂*(1-b)^i
end

"""
    grad_simplex_2D(a, b, id, jd)

Evalute the partial derivative w.r.t. (r,s) of 2D "Legendre" polynomial with
index, or order, (id,jd) at (a,b)

# Examples
```jldoctest
"""

function grad_simplex_2D(a, b, id, jd)
    fa = jacobiP(a, 0, 0, id)
    gb = jacobiP(b, 2*id+1, 0, jd)
    dfa = grad_jacobiP(a, 0, 0, id)
    dgb = grad_jacobiP(b, 2*id+1, 0, jd)

    dmodedr = @. dfa*gb
    if id > 0
        dmodedr = dmodedr.*((0.5*(1 .- b)).^(id-1))
    end

    dmodeds = @. dfa*(gb*(0.5*(1+a)))
    if id > 0
        dmodeds = dmodeds.*((0.5*(1 .- b)).^(id-1))
    end

    tmp = @. dgb*((0.5*(1-b))^id)
    if id > 0
        tmp = tmp-0.5*id*gb.*((0.5*(1 .- b)).^(id-1))
    end
    dmodeds = dmodeds+fa.*tmp

    dmodedr = 2^(id+0.5)*dmodedr
    dmodeds = 2^(id+0.5)*dmodeds
    return dmodedr, dmodeds
end


"""
    rstoab(r, s)

Transfer from reference element coordinate (r,s) to polynomial basis evaluation
coordinate (a,b)

# Examples
```jldoctest
"""
function rstoab(r, s)
    a = zeros(length(r), 1)
    for n = 1:length(r)
        if s[n] != 1
            a[n] = 2*(1+r[n])/(1-s[n])-1
        else
            a[n] = -1
        end
    end
    return a, s
end


"""
    xytors(x, y)

Transfer from equilateral triangle coordinate (a,b) to reference element
coordinate (r,s)

# Examples
```jldoctest
"""
function xytors(x, y)
    L₁ = @. (sqrt(3.0)y+1.0)/3.0
    L₂ = @. (-3.0x-sqrt(3.0)y+2.0)/6.0
    L₃ = @. (3.0x-sqrt(3.0)y+2.0)/6.0
    r = -L₂+L₃-L₁
    s = -L₂-L₃+L₁
    return r, s
end

"""
    function warp_factor(N, rout)

Compute warp factor for triangular interp nodes.

# Examples
```jldoctest
"""

function warp_factor(N, rout)
    LGLr,w1D = gauss_lobatto_quad(0, 0, N)
    # LGLr = transpose(LGLr[:])
    req = LinRange(-1,1,N+1)
    Veq = vandermonde_1D(N, req)
    Nr = length(rout)
    Pmat = zeros(N+1, Nr)
    for i = 1:N+1
        Pmat[i,:] = transpose(jacobiP(rout, 0, 0, i-1))
    end
    Lmat = transpose(Veq)\Pmat
    # warp = transpose(Lmat)*transpose(LGLr.-req)
    warp = transpose(Lmat)*(LGLr.-req[:])
    zerof = @. (abs(rout) < 1.0 - 1.0e-10)
    sf = @. 1.0 - (zerof*rout)^2
    warp = @. warp/sf + warp*(zerof-1)
    return warp
end


"""
    wb_nodes_tri(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function wb_nodes_tri(N)

    Np = convert(Int,(N+1)*(N+2)/2)

    alpopt = [0.0 0.0 1.4152 0.1001 0.2751 0.98 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.577 1.6223 1.6258]
    α = N < 16 ? alpopt[N] : 5/3

    L1 = zeros(Np, 1)
    L2 = zeros(Np, 1)
    L3 = zeros(Np, 1)
    sk = 1
    for n = 1:N+1
        for m = 1:N+2-n
            L1[sk] = (n-1)/N
            L3[sk] = (m-1)/N
            sk += 1;
        end
    end
    L2 = @. 1.0-L1-L3
    x = -L2+L3
    y = (-L2-L3+2*L1)/sqrt(3.0)
    blend1 = 4.0*L2.*L3
    blend2 = 4.0*L1.*L3
    blend3 = 4.0*L1.*L2
    warpf1 = warp_factor(N, L3-L2)
    warpf2 = warp_factor(N, L1-L3)
    warpf3 = warp_factor(N, L2-L1)
    warp1 = @. blend1*warpf1*(1+(α*L1).^2)
    warp2 = @. blend2*warpf2*(1+(α*L2).^2)
    warp3 = @. blend3*warpf3*(1+(α*L3).^2)

    x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3
    y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3

    r, s = xytors(x, y)
    return r[:], s[:]
end

"""
    equi_nodes_tri(N)

Compute optimized interpolation nodes using blend & warp method on equilateral
triangles for polynomial of order N, with Np points

# Examples
```jldoctest
"""

function equi_nodes_tri(N)
    Np = convert(Int,(N+1)*(N+2)/2)

    r = zeros(Np,1)
    s = zeros(Np,1)

    r1D = LinRange(-1,1,N+1)
    sk = 1
    for i = 0:N
        for j = 0:N-i
            r[sk] = r1D[i+1]
            s[sk] = r1D[j+1]
            sk += 1
        end
    end

    return r[:], s[:]
end


"""
    quad_nodes_tri(N)

Returns quadrature nodes which exactly integrate degree N polynomials

# Examples
```jldoctest
"""

function quad_nodes_tri(N)

    if N<28
        rsw = readdlm(string("./src/QuadratureData/quad_nodes_tri_N", N, ".txt"),' ', Float64, '\n')
        r = rsw[:,1]
        s = rsw[:,2]
        w = rsw[:,3]
    else
        cubNA = convert(Int,ceil( (N+1)/2))
        cubNB = convert(Int,ceil( (N+1)/2))
        cubA,cubWA = gauss_quad(0,0, cubNA-1)
        cubB,cubWB = gauss_quad(1,0, cubNB-1)

        cubA = ones(cubNB,1)*cubA'
        cubB = cubB*ones(1,cubNA)

        r = @. 0.5*(1+cubA)*(1-cubB)-1
        s = cubB
        w = 0.5*cubWB*(cubWA')
    end

    return r[:], s[:], w[:]
end
