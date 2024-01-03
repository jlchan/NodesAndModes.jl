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

    @test NodesAndModes.face_vertices(Line()) == (1, 2)
end

area(elem::Quad) = 4.0
area(elem::Tri) = 2.0

@testset "2D $elem basis tests" for elem in (Tri(), Quad())

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

    # elem is either a tri or quad type
    num_faces = (elem == Tri()) ? 3 : 4
    @test length(NodesAndModes.face_vertices(elem)) == num_faces

    if elem == Tri()
        a, b = NodesAndModes.abtors(Tri(), r, s)
        r2, s2 = NodesAndModes.rstoab(Tri(), a, b)
        r3, s3 = NodesAndModes.rstoab(a, b)
        @test r ≈ r2 ≈ r3
        @test s ≈ s2 ≈ s3
    end

    # check for Kronecker structure
    if elem == Quad()
        r, s = nodes(elem, N)
        V = vandermonde(elem, N, r, s)
        V1D = vandermonde(Line(), N, nodes(Line(), N))
        @test V ≈ kron(V1D, V1D)

        rq, sq = quad_nodes(elem, N)
        Vq = vandermonde(elem, N, rq, sq)
        Vq1D = vandermonde(Line(), N, quad_nodes(Line(), N)[1])
        @test Vq ≈ kron(Vq1D, Vq1D)

        @test NodesAndModes.face_basis(Quad(), N, r, s) ≈ NodesAndModes.edge_basis(Quad(), N, r, s)

        # check invertibility of face basis matrix
        Fmask = hcat(NodesAndModes.find_face_nodes(Quad(), r, s)...)
        @test cond(NodesAndModes.face_basis(Quad(), N, r[Fmask], s[Fmask])) < 1e3
    end
end

volume(elem::Hex) = 8
volume(elem::Wedge) = 4
volume(elem::Pyr) = 8/3
volume(elem::Tet) = 4/3

num_faces(::Tet) = 4
num_faces(::Wedge) = 5
num_faces(::Pyr) = 5
num_faces(::Hex) = 6

@testset "3D $elem basis tests" for elem in (Hex(), Wedge(), Pyr(), Tet())

    tol = 5e2 * eps()

    N = 3
    rq, sq, tq, wq = quad_nodes(elem,N)
    @test sum(wq) ≈ volume(elem)
    if elem==Hex()
        @test abs(sum(rq .* wq)) + abs(sum(sq .* wq)) + abs(sum(tq .* wq)) < tol
    end

    Vq = vandermonde(elem, N, rq, sq, tq)
    @test Vq' * diagm(wq) * Vq ≈ I

    @test length(NodesAndModes.face_vertices(elem)) == num_faces(elem)

    if elem != Pyr()
        r, s, t = nodes(elem, N)
        V = vandermonde(elem, N, r, s, t)
        Dr, Ds, Dt = (A->A / V).(grad_vandermonde(elem, N, r, s, t))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) + norm(sum(Dt, dims=2)) < tol
        @test norm(Dr * s) + norm(Dr * t) < tol
        @test norm(Ds * r) + norm(Ds * t) < tol
        @test norm(Dt * r) + norm(Dt * s) < tol
        @test Dr * r ≈ ones(length(r))
        @test Ds * s ≈ ones(length(s))
        @test Dt * t ≈ ones(length(t))

        r, s, t = equi_nodes(elem, N)
        V = vandermonde(elem, N, r, s, t)
        Dr, Ds, Dt = (A->A / V).(grad_vandermonde(elem, N, r, s, t))
        @test norm(sum(Dr, dims=2)) + norm(sum(Ds, dims=2)) + norm(sum(Dt, dims=2)) < tol
        @test norm(Dr * s) + norm(Dr * t) + norm(Ds * r) + norm(Ds * t) + norm(Dt * r) + norm(Dt * s) < tol
        @test Dr * r ≈ one.(r)
        @test Ds * s ≈ one.(s)
        @test Dt * t ≈ one.(t)   
        
    elseif elem == Pyr()
        r, s, t = nodes(elem, N)
        rq, sq, tq, wq = quad_nodes(elem, N)        
        Vq, Vr, Vs, Vt = (A -> A / vandermonde(elem, N, r, s, t)).(basis(elem, N, rq, sq, tq))
        M = Vq' * diagm(wq) * Vq

        # construct weak differentiation matrices instead of nodal differentiation matrices since 
        # the pyramid basis is singular at the node at the tip of the pyramid. 
        Dr, Ds, Dt = (A -> M \ (Vq' * diagm(wq) * A)).((Vr, Vs, Vt))
        @test norm(sum(Dr, dims=2)) < tol
        @test norm(sum(Ds, dims=2)) < tol
        @test norm(sum(Dt, dims=2)) < tol
        @test norm(Dr * s) + norm(Dr * t) < tol
        @test norm(Ds * r) + norm(Ds * t) < tol
        @test norm(Dt * r) + norm(Dt * s) < tol
        @test Dr * r ≈ one.(r)
        @test Ds * s ≈ one.(s)
        @test Dt * t ≈ one.(t)

        a, b, c = NodesAndModes.rsttoabc(Pyr(), rq, sq, tq)
        r2, s2, t2 = NodesAndModes.abctorst(Pyr(), a, b, c)
        @test r2 ≈ rq 
        @test s2 ≈ sq 
        @test t2 ≈ tq 
    end
end

# these can be used for transfinite interpolation
@testset "3D face bases (no interior nodes)" begin

    r, s, t = equi_nodes(Hex(), 3)
    @test size(NodesAndModes.face_basis(Hex(), 2, r, s, t), 2) == 26
    r, s, t = equi_nodes(Tet(), 3)
    @test NodesAndModes.face_basis(Tet(), 2, r, s, t) ≈ NodesAndModes.edge_basis(Tet(), 2, r, s, t)

    @test_throws MethodError NodesAndModes.face_basis(Wedge(), 3, equi_nodes(Wedge(), 4)...)
    @test_throws MethodError NodesAndModes.face_basis(Pyr(), 3, equi_nodes(Pyr(), 4)...)

    N = 3
    for elem in (Hex(), Tet())
        r, s, t = equi_nodes(elem, N)
        V = NodesAndModes.face_basis(elem, N, r, s, t)
        if elem isa Hex
            @test size(V, 2) == 56
        elseif elem isa Tet
            @test size(V, 2) == 20
        end
        @test minimum(svdvals(V)) > 0 # invertibility
        @test norm(V * (V \ (1 .+ r + s + t)) - (1 .+ r + s + t)) < 100 * eps() # polynomial recovery
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

# test high order symmetric quadrature rules
@testset "Simplicial symmetric quadrature" begin
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

# test high order Stroud quadrature
@testset "Simplicial Stroud quadrature" begin
    tol = 1e3*eps()    

    N = 26
    rq, sq, wq = quad_nodes(Tri(), N)
    @test sum(wq) ≈ 2
    @test sum(rq .* wq) ≈ -2/3
    @test sum(sq .* wq) ≈ -2/3

    N = 11
    rq, sq, tq, wq = quad_nodes(Tet(), N)
    @test sum(wq) ≈ 4/3
    @test sum(rq .* wq) ≈ -2/3
    @test sum(sq .* wq) ≈ -2/3
    @test sum(tq .* wq) ≈ -2/3

    @test_throws ArgumentError jaskowiec_sukumar_quad_nodes(Tet(), 21)
end

# apparently inference got better after 1.5 but breaks in nightly
if VERSION >= v"1.6" && VERSION <= v"1.9" 
    @testset "Inferrability tests" begin
        for elementType in (Line(), Tri(), Quad(), Tet(), Pyr(), Wedge(), Hex())
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
end