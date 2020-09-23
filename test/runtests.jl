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

        V = Line.vandermonde(N,r)
        @test V'*diagm(w)*V ≈ I

        Vr = Line.grad_vandermonde(N,r)
        Dr = Vr/V
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)

        r,w = gauss_lobatto_quad(0,0,N)
        @test sum(w) ≈ 2
        @test abs(sum(w.*r)) < tol

        # check if Lobatto is exact for 2N-2 polynoms
        V = Line.vandermonde(N-1,r)
        @test V'*diagm(w)*V ≈ I

        Dr = Line.grad_vandermonde(N,r)/Line.vandermonde(N,r)
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)
    end
end

@testset "2D tri tests" begin
    tol = 1e2*eps()

    N = 3
    rq,sq,wq = Tri.quad_nodes(N)
    @test sum(wq)≈2
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3

    Vq = Tri.vandermonde(N,rq,sq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s = Tri.nodes(N)
    V = Tri.vandermonde(N,r,s)
    Dr,Ds = (A->A/V).(Tri.grad_vandermonde(N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    r,s = Tri.equi_nodes(N)
    V = Tri.vandermonde(N,r,s)
    Dr,Ds = (A->A/V).(Tri.grad_vandermonde(N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    # test duffy quadrature too
    N = 14
    rq,sq,wq = Tri.quad_nodes(N)
    @test sum(wq)≈2
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3
end

@testset "2D quad tests" begin
    tol = 1e2*eps()

    N = 3
    rq,sq,wq = Quad.quad_nodes(N)
    @test sum(wq) ≈ 4
    @test abs(sum(rq.*wq)) < tol
    @test abs(sum(sq.*wq)) < tol

    Vq = Quad.vandermonde(N,rq,sq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s = Quad.nodes(N)
    V = Quad.vandermonde(N,r,s)
    Dr,Ds = (A->A/V).(Quad.grad_vandermonde(N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))

    r,s = Quad.equi_nodes(N)
    V = Quad.vandermonde(N,r,s)
    Dr,Ds = (A->A/V).(Quad.grad_vandermonde(N,r,s))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
    @test norm(Dr*s)+norm(Ds*r) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
end


@testset "3D hex tests" begin
    tol = 5e2*eps()

    N = 3
    rq,sq,tq,wq = Hex.quad_nodes(N)
    @test sum(wq) ≈ 8
    @test abs(sum(rq.*wq))+abs(sum(sq.*wq))+abs(sum(tq.*wq)) < tol

    Vq = Hex.vandermonde(N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = Hex.nodes(N)
    V = Hex.vandermonde(N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(Hex.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    r,s,t = Hex.equi_nodes(N)
    V = Hex.vandermonde(N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(Hex.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))
end


@testset "3D pyr tests" begin
    tol = 5e2*eps()

    N = 3
    rq,sq,tq,wq = Pyr.quad_nodes(N)
    @test sum(wq) ≈ 8/3

    Vq = Pyr.vandermonde(N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = Pyr.nodes(N)
    V = Pyr.vandermonde(N,r,s,t)
    rq,sq,tq,wq = Pyr.quad_nodes(N)
    Vq,Vr,Vs,Vt = (A->A/V).(Pyr.basis(N,rq,sq,tq))
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
    tol = 5e2*eps()

    N = 3
    rq,sq,tq,wq = Wedge.quad_nodes(N)
    @test sum(wq) ≈ 4
    @test abs(sum(tq.*wq)) < tol

    Vq = Wedge.vandermonde(N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    r,s,t = Wedge.nodes(N)
    V = Wedge.vandermonde(N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(Wedge.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    r,s,t = Wedge.equi_nodes(N)
    V = Wedge.vandermonde(N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(Wedge.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))
end

@testset "3D tet tests" begin
    tol = 5e2*eps()

    N = 3
    rq,sq,tq,wq = Tet.quad_nodes(N)
    @test sum(wq) ≈ 4/3

    Vq = Tet.vandermonde(N,rq,sq,tq)
    @test Vq'*diagm(wq)*Vq ≈ I

    # r,s,t = Tet.nodes(N)
    # V = Tet.vandermonde(N,r,s,t)
    # Dr,Ds,Dt = (A->A/V).(Tet.grad_vandermonde(N,r,s,t))
    # @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    # @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    # @test Dr*r ≈ ones(length(r))
    # @test Ds*s ≈ ones(length(s))
    # @test Dt*t ≈ ones(length(t))

    r,s,t = Tet.equi_nodes(N)
    V = Tet.vandermonde(N,r,s,t)
    Dr,Ds,Dt = (A->A/V).(Tet.grad_vandermonde(N,r,s,t))
    @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
    @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
    @test Dr*r ≈ ones(length(r))
    @test Ds*s ≈ ones(length(s))
    @test Dt*t ≈ ones(length(t))

    N = 8
    rq,sq,tq,wq = Tet.quad_nodes(N)
    @test sum(wq)≈4/3
    @test sum(rq.*wq)≈ -2/3
    @test sum(sq.*wq)≈ -2/3
    @test sum(tq.*wq)≈ -2/3
end
