using Test
using MultiAffine

import ManifoldGroupUtils: algebra, rand_lie

import ManifoldsBase as MB
using Manifolds
import Random: rand!, AbstractRNG

import Random

rng = Random.default_rng()




@testset "general $G" for G in [
    MultiDisplacement(3, 2),
    MultiDisplacement(2),
    MultiAffineGroup(Unitary(2), 2),
]
    # the following seed is necessary,
    # for some random cases either of both of these can happen
    # 1. exp_lie is not the inverse of log_lie (due to a bug in julia's matrix log)
    # 2. adjoint action is not exactly an algebra morphism (because of rounding errors)
    Random.seed!(rng, 3)
    n = 3
    pts = [rand(rng, G) for i in 1:n]
    vels = [rand(rng, algebra(G)) for i in 1:n]
    Manifolds.test_group(G, pts, vels, vels,
        test_exp_lie_log=!isa(G, MultiAffineGroup{<:Unitary}),
        test_lie_bracket=true,
        test_adjoint_action=true,
        test_diff=true,
    )
end







@testset "exp (ad_ξ) = Ad_exp(ξ)" for G in [
    MultiDisplacement(3,2),
    MultiDisplacement(2),
    MultiAffineGroup(Unitary(3), 3), # broken: exp(ad_{ξ}) cannot be computed
]
    vel = rand_lie(rng, G)
    tvel = rand_lie(rng, G)
    @test GT.check_exp_ad(G, vel, tvel) broken=G isa MultiAffineGroup{<:Unitary}
end

_adjoint_action(G::MultiAffineGroup, p, X) = begin
    tmp = allocate_result(G, adjoint_action, X)
    return _adjoint_action!(G, tmp, p, X)
end

_adjoint_action!(G::MultiAffineGroup, tmp, p, X) = begin
    mat = affine_matrix(G, p)
    matinv = affine_matrix(G, inv(G, p))
    res = mat * screw_matrix(G, X) * matinv
    map(copyto!, submanifold_components(G, tmp), submanifold_components(G, res))
    return tmp
end

@testset "Compare Adjoint Implementations" for G in [
    MultiDisplacement(3, 2),
    MultiDisplacement(2),
    MultiAffineGroup(Unitary(3), 2),
]
    χ = rand(rng, G)
    ξ = rand_lie(rng, G)
    expected = _adjoint_action(G, χ, ξ)
    computed = adjoint_action(G, χ, ξ)
    @test isapprox(algebra(G), expected, computed)
end


@testset "Diff $G" for G in [
    MultiDisplacement(3, 2),
    SpecialEuclidean(3),
    SpecialOrthogonal(3),
    ]
    @testset "$G $side" for side in [LeftSide(), RightSide()]
        ξ = rand(rng, algebra(G))
        @test GT.check_apply_diff_group_at_id(G, ξ, side) broken=(G==SpecialEuclidean(3))&&(side==RightSide())
    end
end

@testset "Test types" begin
    @testset "MultiDisplacement($dim,$size) creates proper type" for dim in [4]
        for size in [2]
            dim = 4
            size = 2
            GM = MultiDisplacement(dim, size)
            @test isa(GM, MultiDisplacement{dim,size})
            @test !isa(GM, MultiDisplacement{dim,size + 5})
            @test !isa(GM, MultiDisplacement{dim + 1,size})
        end
    end
    @testset "MultiAffineGroup(G, $size) creates proper type" for dim in [4]
        for size in [3]
            GA = MultiAffineGroup(Orthogonal(dim), size)
            @test isa(GA, MultiAffineGroup{typeof(Orthogonal(dim)),dim,size,ℝ})
            @test !isa(GA, MultiAffineGroup{typeof(Orthogonal(dim + 1)),dim + 1,size,ℝ})
            @test isa(MultiAffineGroup(Unitary(dim), size), MultiAffineGroup)
        end
    end
end


check_from_normal_grp(G::MultiAffineGroup, ts) = begin
    χ1 = from_normal_grp(G, eachcol(ts)...)
    χ2 = from_normal_grp(G, ts)
    return isapprox(G, χ1, χ2)
end

check_from_normal_alg(G::MultiAffineGroup, ts) = begin
    ξ1 = from_normal_alg(G, eachcol(ts)...)
    ξ2 = from_normal_alg(G, ts)
    return isapprox(algebra(G), ξ1, ξ2)
end

"""
    rand_trans(rng, G::MultiAffineGroup)

Random translation part of the group `G`.
"""
rand_trans(rng, G::MultiAffineGroup{TH, dim, size}) where {TH, dim, size} = randn(rng, dim, size)

@testset "from/to $G" for G in [MultiDisplacement(3,2)]
    @testset "grp" begin
        @test check_from_normal_grp(G, rand_trans(rng, G))
    end
    @testset "alg" begin
        @test check_from_normal_alg(G, rand_trans(rng, G))
    end
end


@testset "Test $G" for G in [
    MultiDisplacement(3,2),
    MultiDisplacement(2),
    MultiAffineGroup(Unitary(3), 2),
    ]
    @test GT.check_grp_rep_Identity(G, affine_matrix)
    vel = rand_lie(rng, G)
    pt = rand(rng, G)
    @test GT.check_exp_lie_point(G, vel)
    @test GT.check_adjoint_action_in_alg(G, pt, vel)
    @testset "zero_element" begin
        @test GT.check_zero_Identity(G)
    end
    @testset "Inverse & Matrix" begin
        χ = rand(rng, G)
        @test GT.check_inv_rep(G, affine_matrix, χ)
    end
    @testset "Adjoint action & matrix" begin
        χ = rand(rng, G)
        ξ = rand_lie(rng, G)
        @test GT.check_adjoint_action(G, affine_matrix, screw_matrix, χ, ξ)
    end
    @testset "Lie Bracket & matrix" begin
        v1, v2 = [rand_lie(rng, G) for i in 1:2]
        @test GT.check_alg_rep(G, screw_matrix, v1, v2)
    end
    @testset "Composition & matrix" begin
        p1, p2 = [rand(rng, G) for i in 1:2]
        @test GT.check_grp_rep_compose(G, affine_matrix, p1, p2)
    end
end


@testset for G in [MultiDisplacement(3,2)]
    χ1, χ2 = [rand(rng, G) for i in 1:2]
    ξ1, ξ2 = [rand_lie(rng, G) for i in 1:2]
    v1 = translate_from_id(G, χ1, ξ1, LeftSide())
    @test GT.check_exp_invariant(G, exp, χ1, v1, χ2)
    @test GT.check_exp_log(G, exp, log, χ1, χ2)
    @test GT.check_log_exp(G, log, exp, χ1, v1)
    v = similar(v1)
    @test GT.check_log_log_(G, log, log!, v, χ1, χ2)
    χ = similar(χ1)
    @test GT.check_exp_exp_(G, exp, exp!, χ, χ1, v1)
end



include("multiaffine/apply_diff_group.jl")
include("multiaffine/inv_diff.jl")
