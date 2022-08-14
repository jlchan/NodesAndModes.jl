function basis(elem::Hex, N, r, s, t)
    Np = convert(Int, (N+1)^3)
    V, Vr, Vs, Vt = ntuple(x->zeros(length(r), Np), 4)
    tmp_array = zeros(N+1)
    for ii in eachindex(r, s, t)                                                
        sk = 1
        for k in 0:N
            P_k = jacobiP(t[ii], 0, 0, k)
            dP_k = grad_jacobiP(t[ii], 0, 0, k)
            for i in 0:N
                P_i = jacobiP(r[ii], 0, 0, i)
                dP_i = grad_jacobiP(r[ii], 0, 0, i)
                for j in 0:N
                    P_j = jacobiP(s[ii], 0, 0, j)
                    dP_j = grad_jacobiP(s[ii], 0, 0, j)
                    
                    V[ii, sk]  = P_i * P_j * P_k
                    Vr[ii, sk] = dP_i * P_j * P_k
                    Vs[ii, sk] = P_i * dP_j * P_k
                    Vt[ii, sk] = P_i * P_j * dP_k
                   
                    sk += 1
                end                
            end
        end
    end

    return V, Vr, Vs, Vt
end


# ===================================================

function nodes(elem::Hex, N)
    r1D, _ = gauss_lobatto_quad(0,0,N)
    return vec.(meshgrid(r1D,r1D,r1D))
end


function equi_nodes(elem::Hex, N)
    r1D = LinRange(-1,1,N+1)
    return vec.(meshgrid(r1D,r1D,r1D))
end


function quad_nodes(elem::Hex, N)
    r1D,w1D = gauss_quad(0,0,N)
    r,s,t = vec.(meshgrid(r1D,r1D,r1D))
    wr,ws,wt = vec.(meshgrid(w1D,w1D,w1D))
    w = @. wr*ws*wt
    return r,s,t,w
end
