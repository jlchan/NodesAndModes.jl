function basis(elem::Hex, N, r, s, t)
    Np = convert(Int,(N+1)^3)
    sk = 1
    V, Vr, Vs, Vt = ntuple(x->zeros(length(r), Np),4)
    for i=0:N
        for j=0:N
            for k=0:N
                V[:,sk]  .= jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t, 0, 0, k)
                Vr[:,sk] .= grad_jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vs[:,sk] .= jacobiP(r, 0, 0, i).*grad_jacobiP(s, 0, 0, j).*jacobiP(t,0,0,k)
                Vt[:,sk] .= jacobiP(r, 0, 0, i).*jacobiP(s, 0, 0, j).*grad_jacobiP(t,0,0,k)
                sk += 1
            end
        end
    end

    return V, Vr, Vs, Vt
end

# ===================================================

function nodes(elem::Hex, N)
    r1D,w1D = gauss_lobatto_quad(0,0,N)
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
