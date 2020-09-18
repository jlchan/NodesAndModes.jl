function build_edge_basis(N, basis_1D, vertices,
                          edges, vertex_functions)
    
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
