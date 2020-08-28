using NodesAndModes
using Test
using LinearAlgebra

@testset "1D Basis" begin
    tol = 1e2*eps()

    N = 2
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

    V = vandermonde_1D(N-1,r)
    @test V'*diagm(w)*V ≈ I

    Dr = grad_vandermonde_1D(N,r)/vandermonde_1D(N,r)    
    @test norm(sum(Dr,dims=2)) < tol
    @test Dr*r ≈ ones(N+1)
end
