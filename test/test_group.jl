using MultiAffine

import ManifoldGroupUtils: algebra, rand_lie, translate_from_id

import ManifoldsBase as MB
using Manifolds
import Random: rand!, AbstractRNG

import Random

rng = Random.default_rng()

group_list = [
    MultiDisplacementGroup(3, 2),
    MultiDisplacementGroup(2),
    MultiAffineGroup(Unitary(2), 2),
    MultiAffineGroup(Unitary(3), 2),
]


@testset "general $G" for (i, G) in enumerate(group_list)
    N = length(group_list)
    @info "$i / $N: $G"
    # the following seed is necessary,
    # for some random cases either of both of these can happen
    # 1. exp_lie is not the inverse of log_lie (due to a bug in julia's matrix log)
    # 2. adjoint action is not exactly an algebra morphism (because of rounding errors)
    Random.seed!(rng, 0)
    n = 3
    pts = map(ξ -> exp_lie(G, ξ), [rand_lie(rng, G) for i in 1:n])
    vels = [rand_lie(rng, G) for i in 1:n]
    Manifolds.test_group(
        G, pts, vels, vels,
        test_exp_lie_log=!isa(G, MultiAffineGroup{<:Unitary}),
        test_lie_bracket=false, # TODO: remove
        test_adjoint_action=true,
        test_diff=true,
    )
    @info "–"^16
end







@testset "exp (ad_ξ) = Ad_exp(ξ)" for G in [
    MultiDisplacementGroup(3,2),
    MultiDisplacementGroup(2),
    MultiAffineGroup(Unitary(3), 3), # broken: exp(ad_{ξ}) cannot be computed
    ]
    vel = rand_lie(rng, G)
    tvel = rand_lie(rng, G)
    Test.@test GT.check_exp_ad(G, vel, tvel) broken=G isa MultiAffineGroup{<:Unitary}
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

@testset "Compare Adjoint Implementations" for G in group_list
    χ = rand(rng, G)
    ξ = rand_lie(rng, G)
    expected = _adjoint_action(G, χ, ξ)
    computed = adjoint_action(G, χ, ξ)
    @test isapprox(algebra(G), expected, computed)
end


@testset "Diff $G" for G in group_list
    @testset "$G $side" for side in [LeftSide(), RightSide()]
        ξ = rand(rng, algebra(G))
        Test.@test GT.check_apply_diff_group_at_id(G, ξ, side) broken=(G==SpecialEuclidean(3))&&(side==RightSide())&&(MultiAffine._MANIFOLDS_VERSION < v"0.10")
    end
end

@testset "Test types" begin
    @testset "MultiDisplacementGroup($dim,$size) creates proper type" for dim in [4]
        for size in [2]
            dim = 4
            size = 2
            GM = MultiDisplacementGroup(dim, size)
            @test isa(GM, MultiDisplacementGroup{dim,size})
            @test !isa(GM, MultiDisplacementGroup{dim,size + 5})
            @test !isa(GM, MultiDisplacementGroup{dim + 1,size})
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
rand_trans(rng, G::MultiAffineGroup{<:Any,dim,size}) where {dim,size} = randn(rng, dim, size)

@testset "from/to $G" for G in [MultiDisplacementGroup(3, 2)]
    @testset "grp" begin
        @test check_from_normal_grp(G, rand_trans(rng, G))
    end
    @testset "alg" begin
        @test check_from_normal_alg(G, rand_trans(rng, G))
    end
end

check_indices(G::MultiAffineGroup, v, proj, ind) = begin
    w = proj(G, v)
    idx = eachindex(v)
    return vec(w) == vec(v)[ind(G, idx)]
end

check_proj_point(G, subman, proj, χ) = is_point(subman, proj(G, χ))

@testset "Proj/Indices $G" for G in
    [MultiDisplacementGroup(3,2)]
    χ = rand(rng, G)
    @test check_indices(G, χ, MultiAffine.to_normal, MultiAffine.normal_indices)
    @test check_indices(G, χ, MultiAffine.to_factor, MultiAffine.factor_indices)
    @test check_proj_point(G, factor_group(G), to_factor_grp, χ)
    @test check_proj_point(G, normal_group(G), to_normal_grp, χ)
    ξ = rand_lie(rng, G)
    @test check_proj_point(G, algebra(normal_group(G)), to_normal_alg, ξ)
    @test check_proj_point(G, algebra(factor_group(G)), to_factor_alg, ξ)
    @testset for pos in 0:1
        ind_pos(G, idx) = MultiAffine.normal_indices_at(G, idx, pos)
        proj_pos(G, x) = MultiAffine.to_normal(G, x)[:,pos+1]
        @test check_indices(G, χ, proj_pos, ind_pos)
    end
end


@testset "Test $G" for G in group_list
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


run_test(G, name, args) = begin
    @testset "$name" begin
        @test getfield(GT, Symbol("check_"*name))(G, args...)
    end
end

grp_rep(G::MultiAffineGroup, x) = affine_matrix(G, x)
# grp_rep(G, ::Identity) = grp_rep(G, identity_element(G))
alg_rep(G::MultiAffineGroup, x) = screw_matrix(G, x)

@testset "Group Testing.jl $G" for G in group_list
    χ1, χ2 = map(ξ -> exp_lie(G, ξ), [rand_lie(rng, G) for i in 1:2])
    # χ1, χ2 = [rand(rng, G) for i in 1:2]
    ξ1, ξ2 = [rand_lie(rng, G) for i in 1:2]
    v1 = translate_from_id(G, χ1, ξ1, LeftSide())
    v = similar(v1)
    χ = similar(χ1)
    tests = [
        ("exp_lie_point", [ξ1]),
        ("adjoint_action_in_alg", [χ1, ξ1]),
        ("adjoint_action_lie_bracket", [χ1, ξ1, ξ2]),
        ("grp_rep_Identity", [grp_rep]),
        ("grp_rep_compose", [grp_rep, χ1, χ2]),
        ("alg_rep", [alg_rep, ξ1, ξ2]),
        ("adjoint_action", [grp_rep, alg_rep, χ1, ξ1]),
        ("inv_rep", [grp_rep, χ1]),
        ("zero_Identity", []),
        ("exp_log", [exp, log, χ1, χ2]),
        ("log_exp", [log, exp, χ1, v1]),
        ("log_log_", [log, log!, v, χ1, χ2]),
        ("exp_exp_", [exp, exp!, χ, χ1, v1]),
        ("exp_lie_ad", [χ1, ξ1]),
    ]
    for (name, args) in tests
        # TODO: remove
        if (G isa MultiDisplacementGroup{2}) && (name == "adjoint_action_lie_bracket")
            continue
        end
        run_test(G, name, args)
    end

    if G isa AbstractDecoratorManifold{ℂ}
        @test_throws MethodError GT.check_exp_ad(G, ξ1, ξ2)
    else
        run_test(G, "exp_ad", [ξ1, ξ2])
    end
    @testset "$side" for side in [LeftSide(), RightSide()]
        side_tests = [
            ("apply_diff_group_at_id", [ξ1, side, Identity]),
            ("apply_diff_group_at_id", [ξ1, side, identity_element]),
            ("inv_diff", [χ1, ξ1, side]),
        ]
        for (name, args) in side_tests
            run_test(G, name, args)
        end
    end
    @testset "$conv" for conv in [(LeftAction(), LeftSide()), (RightAction(), RightSide())]
        run_test(G, "exp_invariant", [exp, χ1, v1, χ2, conv])
    end
end


include("multiaffine/apply_diff_group.jl")
include("multiaffine/inv_diff.jl")
