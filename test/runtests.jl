using NodesAndModes
using Test
using LinearAlgebra

@testset "1D tests" begin
    tol = 1e2*eps()

    N = 0
    r,w = gauss_lobatto_quad(0,0,N)
    @test sum(w) ≈ 2

    for N = 1:2
        r,w = gauss_quad(0,0,N)
        @test sum(w) ≈ 2
        @test abs(sum(w.*r)) < tol

        V = vandermonde(Line(),N,r)
        @test V'*diagm(w)*V ≈ I

        Vr = grad_vandermonde(Line(),N,r)
        Dr = Vr/V
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)

        r,w = gauss_lobatto_quad(0,0,N)
        @test sum(w) ≈ 2
        @test abs(sum(w.*r)) < tol

        # check if Lobatto is exact for 2N-2 polynoms
        V = vandermonde(Line(),N-1,r)
        @test V'*diagm(w)*V ≈ I

        Dr = grad_vandermonde(Line(),N,r)/vandermonde(Line(),N,r)
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)
    end
end

@testset "2D tri tests" begin
    tol = 1e2*eps()

    N = 3
    rq,sq,wq = quad_nodes(Tri(),N)
    @test sum(wq)≈2
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3

    Vq = vandermonde(Tri(),N,rq,sq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s = nodes(Tri(),N)
    V = vandermonde(Tri(),N,r,s)
    Dr,Ds = (A->A/V).(grad_vandermonde(Tri(),N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    r,s = equi_nodes(Tri(),N)
    V = vandermonde(Tri(),N,r,s)
    Dr,Ds = (A->A/V).(grad_vandermonde(Tri(),N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    # test duffy quadrature too
    N = 14
    rq,sq,wq = quad_nodes(Tri(),N)
    @test sum(wq)≈2
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3
end

@testset "2D quad tests" begin
    tol = 1e2*eps()

    N = 3
    rq,sq,wq = quad_nodes(Quad(),N)
    @test sum(wq) ≈ 4
    @test abs(sum(rq.*wq)) < tol
    @test abs(sum(sq.*wq)) < tol

    Vq = vandermonde(Quad(),N,rq,sq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s = nodes(Quad(),N)
    V = vandermonde(Quad(),N,r,s)
    Dr,Ds = (A->A/V).(grad_vandermonde(Quad(),N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    r,s = equi_nodes(Quad(),N)
    V = vandermonde(Quad(),N,r,s)
    Dr,Ds = (A->A/V).(grad_vandermonde(Quad(),N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
end

@testset "3D hex tests" begin
    tol = 5e2*eps()

    N = 3
    rq,sq,tq,wq = quad_nodes(Hex(),N)
    @test sum(wq) ≈ 8
    @test abs(sum(rq.*wq))+abs(sum(sq.*wq))+abs(sum(tq.*wq)) < tol

    Vq = vandermonde(Hex(),N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = nodes(Hex(),N)
    V = vandermonde(Hex(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Hex(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    r,s,t = equi_nodes(Hex(),N)
    V = vandermonde(Hex(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Hex(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))
end

@testset "3D pyr tests" begin
    tol = 1e3*eps()

    N = 3
    rq,sq,tq,wq = quad_nodes(Pyr(),N)
    @test sum(wq) ≈ 8/3

    Vq = vandermonde(Pyr(),N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = nodes(Pyr(),N)
    V = vandermonde(Pyr(),N,r,s,t)
    rq,sq,tq,wq = quad_nodes(Pyr(),N)
    Vq,Vr,Vs,Vt = (A->A/V).(basis(Pyr(),N,rq,sq,tq))
    M = Vq'*diagm(wq)*Vq
    Dr,Ds,Dt = (A->M\(Vq'*diagm(wq)*A)).((Vr,Vs,Vt))
    # Dr,Ds,Dt = (A->A/V).(Pyr.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))
end

@testset "3D wedge tests" begin
    tol = 1e3*eps()

    N = 3
    rq,sq,tq,wq = quad_nodes(Wedge(),N)
    @test sum(wq) ≈ 4
    @test abs(sum(tq.*wq)) < tol

    Vq = vandermonde(Wedge(),N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = nodes(Wedge(),N)
    V = vandermonde(Wedge(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Wedge(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    r,s,t = equi_nodes(Wedge(),N)
    V = vandermonde(Wedge(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Wedge(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))
end

@testset "3D tet tests" begin
    tol = 1e3*eps()

    N = 3
    rq,sq,tq,wq = quad_nodes(Tet(),N)
    @test sum(wq) ≈ 4/3

    Vq = vandermonde(Tet(),N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = nodes(Tet(),N)
    V = vandermonde(Tet(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Tet(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    r,s,t = equi_nodes(Tet(),N)
    V = vandermonde(Tet(),N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(grad_vandermonde(Tet(),N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    # test high order quadrature
    N = 8
    rq,sq,tq,wq = quad_nodes(Tet(),N)
    @test sum(wq)≈4/3
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3
    @test sum(tq.*wq)≈ -2/3
end

@testset "Inferrability for $elementType" for elementType in [Line(), Tri(), Quad(), Tet(), Pyr(), Wedge(), Hex()] 
    N = 1
    if elementType==Line()
        @test (@inferred gauss_quad(0,0,N)) == ([-0.5773502691896257, 0.5773502691896257], [0.9999999999999998, 0.9999999999999998])
        @test (@inferred gauss_lobatto_quad(0,0,N)) == ([-1.0, 1.0], [1.0, 1.0])
    end
    
    @inferred nodes(elementType,N)
    @inferred quad_nodes(elementType,N)
    @inferred equi_nodes(elementType,N)
    if elementType==Line()
        @inferred basis(elementType,N,nodes(elementType,N))
    else
        @inferred basis(elementType,N,nodes(elementType,N)...)
    end
end