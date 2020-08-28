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

        V = vandermonde_1D(N,r)
        @test V'*diagm(w)*V ≈ I

        Vr = grad_vandermonde_1D(N,r)
        Dr = Vr/V
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)

        r,w = gauss_lobatto_quad(0,0,N)
        @test sum(w) ≈ 2
        @test abs(sum(w.*r)) < tol

        # check if Lobatto is exact for 2N-2 polynoms
        V = vandermonde_1D(N-1,r)
        @test V'*diagm(w)*V ≈ I

        Dr = grad_vandermonde_1D(N,r)/vandermonde_1D(N,r)
        @test norm(sum(Dr,dims=2)) < tol
        @test Dr*r ≈ ones(N+1)
    end
end

# @testset "2D tri tests" begin
#     tol = 1e2*eps()
#
#     N = 3
#     rq,sq,wq = Tri.quad_nodes_2D(2*N)
#     @test sum(wq)≈2
#     @test sum(rq.*wq)≈ -2/3
#     @test sum(sq.*wq)≈ -2/3
#
#     Vq = Tri.vandermonde_2D(N,rq,sq)
#     @test Vq'*diagm(wq)*Vq ≈ I
#
#     r,s = Tri.nodes_2D(N)
#     V = Tri.vandermonde_2D(N,r,s)
#     Dr,Ds = (A->A/V).(Tri.grad_vandermonde_2D(N,r,s))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
#     @test norm(Dr*s)+norm(Ds*r) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
#
#     r,s = Tri.equi_nodes_2D(N)
#     V = Tri.vandermonde_2D(N,r,s)
#     Dr,Ds = (A->A/V).(Tri.grad_vandermonde_2D(N,r,s))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
#     @test norm(Dr*s)+norm(Ds*r) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
#
#     # test duffy quadrature too
#     N = 15
#     rq,sq,wq = Tri.quad_nodes_2D(2*N)
#     @test sum(wq)≈2
#     @test sum(rq.*wq)≈ -2/3
#     @test sum(sq.*wq)≈ -2/3
# end
#
# @testset "2D quad tests" begin
#     tol = 1e2*eps()
#
#     N = 3
#     rq,sq,wq = Quad.quad_nodes_2D(2*N)
#     @test sum(wq) ≈ 4
#     @test abs(sum(rq.*wq)) < tol
#     @test abs(sum(sq.*wq)) < tol
#
#     Vq = Quad.vandermonde_2D(N,rq,sq)
#     @test Vq'*diagm(wq)*Vq ≈ I
#
#     r,s = Quad.nodes_2D(N)
#     V = Quad.vandermonde_2D(N,r,s)
#     Dr,Ds = (A->A/V).(Quad.grad_vandermonde_2D(N,r,s))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
#     @test norm(Dr*s)+norm(Ds*r) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
#
#     r,s = Quad.equi_nodes_2D(N)
#     V = Quad.vandermonde_2D(N,r,s)
#     Dr,Ds = (A->A/V).(Quad.grad_vandermonde_2D(N,r,s))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) < tol
#     @test norm(Dr*s)+norm(Ds*r) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
# end
#
#
# @testset "3D hex tests" begin
#     tol = 5e2*eps()
#
#     N = 3
#     rq,sq,tq,wq = Hex.quad_nodes_3D(2*N)
#     @test sum(wq) ≈ 8
#     @test abs(sum(rq.*wq))+abs(sum(sq.*wq))+abs(sum(tq.*wq)) < tol
#
#     Vq = Hex.vandermonde_3D(N,rq,sq,tq)
#     @test Vq'*diagm(wq)*Vq ≈ I
#
#     r,s,t = Hex.nodes_3D(N)
#     V = Hex.vandermonde_3D(N,r,s,t)
#     Dr,Ds,Dt = (A->A/V).(Hex.grad_vandermonde_3D(N,r,s,t))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
#     @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
#     @test Dt*t ≈ ones(length(t))
#
#     r,s,t = Hex.equi_nodes_3D(N)
#     V = Hex.vandermonde_3D(N,r,s,t)
#     Dr,Ds,Dt = (A->A/V).(Hex.grad_vandermonde_3D(N,r,s,t))
#     @test norm(sum(Dr,dims=2)) + norm(sum(Ds,dims=2)) + norm(sum(Dt,dims=2)) < tol
#     @test norm(Dr*s)+norm(Dr*t)+norm(Ds*r)+norm(Ds*t)+norm(Dt*r)+norm(Dt*s) < tol
#     @test Dr*r ≈ ones(length(r))
#     @test Ds*s ≈ ones(length(s))
#     @test Dt*t ≈ ones(length(t))
# end
