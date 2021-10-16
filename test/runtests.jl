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

@testset "2D $elem basis tests" for elem in (Tri(), Quad())
    area(elem::Quad) = 4.0
    area(elem::Tri) = 2.0

    @testset begin
        tol = 1e2*eps()

        N = 3
        rq, sq, wq = quad_nodes(elem, N)
        @test sum(wq) ≈ area(elem) 

        Vq = vandermonde(elem, N, rq, sq)
        @test Vq' * diagm(wq) * Vq ≈ I

        r, s = nodes(elem, N)
        V = vandermonde(elem, N, r, s)
        Dr, Ds = (A->A / V).(grad_vandermonde(elem, N, r, s))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) < tol
        @test norm(Dr * s) + norm(Ds * r) < tol
        @test Dr * r ≈ ones(length(r))
        @test Ds * s ≈ ones(length(s))

        r, s = equi_nodes(elem, N)
        V = vandermonde(elem, N, r, s)
        Dr,Ds = (A->A / V).(grad_vandermonde(elem, N, r, s))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) < tol
        @test norm(Dr * s) + norm(Ds * r) < tol
        @test Dr*r ≈ ones(length(r))
        @test Ds*s ≈ ones(length(s))
    end
end

@testset "3D $elem basis tests" for elem in (Hex(), Wedge(), Pyr(), Tet())

    volume(elem::Hex) = 8
    volume(elem::Wedge) = 4
    volume(elem::Pyr) = 8/3
    volume(elem::Tet) = 4/3

    @testset "3D tests for $elem" begin
        tol = 5e2*eps()

        N = 3
        rq, sq, tq, wq = quad_nodes(elem,N)
        @test sum(wq) ≈ volume(elem)
        if elem==Hex()
            @test abs(sum(rq .* wq)) + abs(sum(sq .* wq)) + abs(sum(tq .* wq)) < tol
        end

        Vq = vandermonde(elem, N, rq, sq, tq)
        @test Vq' * diagm(wq) * Vq ≈ I

        r, s, t = nodes(elem, N)
        V = vandermonde(elem, N, r, s, t)
        Dr, Ds, Dt = (A->A / V).(grad_vandermonde(elem,N, r, s, t))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) + norm(sum(Dt, dims=2)) < tol
        @test norm(Dr * s) + norm(Dr * t) + norm(Ds * r) + norm(Ds * t) + norm(Dt * r) + norm(Dt * s) < tol
        @test Dr*r ≈ ones(length(r))
        @test Ds*s ≈ ones(length(s))
        @test Dt*t ≈ ones(length(t))

        r, s, t = equi_nodes(elem, N)
        V = vandermonde(elem,N, r, s, t)
        Dr, Ds, Dt = (A->A / V).(grad_vandermonde(elem, N, r, s, t))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) + norm(sum(Dt, dims=2)) < tol
        @test norm(Dr * s) + norm(Dr * t) + norm(Ds * r) + norm(Ds * t) + norm(Dt * r) + norm(Dt * s) < tol
        @test Dr * r ≈ ones(length(r))
        @test Ds * s ≈ ones(length(s))
        @test Dt * t ≈ ones(length(t))
    end
end

@testset "Test for Kronecker structure in the Hex basis matrix" begin
    tol = 5e2*eps()
    N = 3
    VDM_1D = vandermonde(Line(), N, nodes(Line(), N))
    VDM_quad_1D = vandermonde(Line(), N, quad_nodes(Line(), N)[1])
    @test abs(norm(vandermonde(Hex(), N, nodes(Hex(), N)...) - kron(VDM_1D, VDM_1D, VDM_1D))) < tol
    @test abs(norm(vandermonde(Hex(), N, quad_nodes(Hex(), N)[1:3]...) - kron(VDM_quad_1D, VDM_quad_1D, VDM_quad_1D))) < tol
end

# test high order quadrature
@testset "Simplicial Stroud quadrature" begin
    tol = 1e3*eps()

    N = 14
    rq, sq, wq = quad_nodes(Tri(), N)
    @test sum(wq) ≈ 2
    @test sum(rq .* wq) ≈ -2/3
    @test sum(sq .* wq) ≈ -2/3

    N = 8
    rq, sq, tq, wq = quad_nodes(Tet(), N)
    @test sum(wq) ≈ 4/3
    @test sum(rq .* wq) ≈ -2/3
    @test sum(sq .* wq) ≈ -2/3
    @test sum(tq .* wq) ≈ -2/3
end

if VERSION >= v"1.6" # apparently inference got better after 1.5?
    @testset "Inferrability for $elementType" for elementType in (Line(), Tri(), Quad(), Tet(), Pyr(), Wedge(), Hex())
        N = 1
        if elementType == Line()
            @test_nowarn (@inferred gauss_quad(0, 0, N)) == ([-0.5773502691896257, 0.5773502691896257], [0.9999999999999998, 0.9999999999999998])
            @test_nowarn (@inferred gauss_lobatto_quad(0, 0, N)) == ([-1.0, 1.0], [1.0, 1.0])
        end
        
        @test_nowarn @inferred nodes(elementType, N)
        @test_nowarn @inferred quad_nodes(elementType, N)
        @test_nowarn @inferred equi_nodes(elementType, N)
        if elementType == Line()
            @test_nowarn @inferred basis(elementType, N, nodes(elementType,N))
        else
            @test_nowarn @inferred basis(elementType, N, nodes(elementType,N)...)
        end
    end
end