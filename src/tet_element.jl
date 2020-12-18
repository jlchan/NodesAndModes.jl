#####
##### 3D modes on tetrahedra
#####

"""
    simplex_3D(a, b, c, i, j, k)

Evaluate 3D "Legendre" basis phi_ijk at (a,b,c) coordinates on the [-1,1] cube
"""
function simplex_3D(a, b, c, i, j, k)
    h1 = jacobiP(a,0,0,i)
    h2 = jacobiP(b,2*i+1,0,j)
    h3 = jacobiP(c,2*(i+j)+2,0,k)
    return @. 2*sqrt(2)*h1*h2*((1-b)^i)*h3*((1-c)^(i+j))
end

"""
    grad_simplex_3D(a, b, c, id, jd, kd)

Evalute the partial derivative w.r.t. (r,s,t) of 3D "Legendre" polynomial with
index, or order, (id,jd,kd) at (a,b,c)
"""

function grad_simplex_3D(a, b, c, id, jd, kd)
    fa = jacobiP(a,0,0,id)
    gb = jacobiP(b,2*id+1,0,jd)
    hc = jacobiP(c,2*(id+jd)+2,0,kd)

    dfa = grad_jacobiP(a,0,0,id)
    dgb = grad_jacobiP(b,2*id+1,0,jd)
    dhc = grad_jacobiP(c,2*(id+jd)+2,0,kd)

    # r-derivative
    V3Dr = @. dfa.*(gb.*hc)
    if id>0
        V3Dr = @. V3Dr*((0.5*(1-b))^(id-1))
    end
    if id+jd>0
        V3Dr = @. V3Dr*((0.5*(1-c))^(id+jd-1))
    end

    # s-derivative
    V3Ds = @. 0.5*(1+a)*V3Dr
    tmp = @. dgb*((0.5*(1-b))^id)
    if id>0
        tmp = @. tmp+(-0.5*id)*(gb*(0.5*(1-b))^(id-1))
    end
    if id+jd>0
        tmp = @. tmp*((0.5*(1-c))^(id+jd-1))
    end
    tmp = @. fa*(tmp*hc)
    V3Ds = @. V3Ds+tmp

    # t-derivative
    V3Dt = @. 0.5*(1+a).*V3Dr+0.5*(1+b)*tmp
    tmp = @. dhc*((0.5*(1-c))^(id+jd))
    if id+jd>0
      tmp = @. tmp-0.5*(id+jd)*(hc*((0.5*(1-c))^(id+jd-1)))
    end
    tmp = @. fa*(gb*tmp)
    tmp = @. tmp*((0.5*(1-b))^id)
    V3Dt = @. V3Dt+tmp

    # normalize
    V3Dr = @. V3Dr*(2^(2*id+jd+1.5))
    V3Ds = @. V3Ds*(2^(2*id+jd+1.5))
    V3Dt = @. V3Dt*(2^(2*id+jd+1.5))
    return V3Dr,V3Ds,V3Dt
end


function equi_nodes(elem::Tet,N)
    Np = (N+1)*(N+2)*(N+3) ÷ 6
    r1D = LinRange(-1,1,N+1)
    r,s,t = ntuple(x->zeros(Np),3)
    sk = 1
    for i = 0:N
        for j = 0:N-i
            for k = 0:N-i-j
                r[sk] = r1D[i+1]
                s[sk] = r1D[j+1]
                t[sk] = r1D[k+1]
                sk += 1
            end
        end
    end

    return r,s,t
end


"""
    quad_nodes_tet(N)

Returns quadrature nodes and weights which exactly integrate degree N polynomials
"""
function quad_nodes_tet(N)

    if N<16
        rstw = readdlm(string(@__DIR__,"/QuadratureData/quad_nodes_tet_N", N, ".txt"),' ', Float64, '\n')
        r = rstw[:,1]
        s = rstw[:,2]
        t = rstw[:,3]
        w = rstw[:,4]
    else
        cubN = convert(Int,ceil((N+1)/2))
        r,s,t,w = stroud_quad_nodes(Tet(),cubN)
    end

    return r,s,t,w
end

function stroud_quad_nodes(elem::Tet,N)
    cubA,cubWA = gauss_quad(0,0,N)
    cubB,cubWB = gauss_quad(1,0,N)
    cubC,cubWC = gauss_quad(2,0,N)

    a,b,c = vec.(meshgrid(cubA,cubB,cubC))
    wa,wb,wc = vec.(meshgrid(cubWA,cubWB,cubWC))

    r = @. .5*(1+a)*.5*(1-b)*(1-c)-1
    s = @. .5*(1+b)*(1-c)-1
    t = c
    w = @. wa*wb*wc
    w = (4/3)*w./sum(w)
    return r,s,t,w
end


quad_nodes(elem::Tet,N) = quad_nodes_tet(2*N)

function basis(elem::Tet,N,r,s,t)
    Np = (N+1)*(N+2)*(N+3)÷6

    V,Vr,Vs,Vt = ntuple(x->zeros(length(r),Np),4)

    a = @. 2*(1+r)/(-s-t)-1
    b = @. 2*(1+s)/(1-t)-1
    c = t

    tol = 1e-14
    ida = @. abs(s+t) < tol
    idb = @. abs(1-t) < tol
    a[ida] .= -1
    b[idb] .= -1

    sk = 1
    for i = 0:N
        for j = 0:N-i
            for k = 0:N-i-j
                V[:,sk] = simplex_3D(a,b,c,i,j,k)
                Vr[:,sk],Vs[:,sk],Vt[:,sk] = grad_simplex_3D(a,b,c,i,j,k)
                sk += 1
            end
        end
    end
    return V,Vr,Vs,Vt
end
