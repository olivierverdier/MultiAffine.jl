using Test
using MultiAffinity

import GroupTools

import ManifoldsBase as MB
using Manifolds
import Random: rand!, AbstractRNG

import Random

rng = Random.default_rng()





@testset "general $G" for G in [
    MultiDisplacement(3, 2),
    MultiDisplacement(2),
    MultiAffine(Unitary(2), 2),
]
    # the following seed is necessary,
    # for some random cases either of both of these can happen
    # 1. exp_lie is not the inverse of log_lie (due to a bug in julia's matrix log)
    # 2. adjoint action is not exactly an algebra morphism (because of rounding errors)
    Random.seed!(rng, 3)
    n = 3
    pts = [rand(rng, G) for i in 1:n]
    vels = [rand(rng, GroupTools.algebra(G)) for i in 1:n]
    Manifolds.test_group(G, pts, vels, vels,
        test_exp_lie_log=!isa(G, MultiAffine{<:Unitary}),
        test_lie_bracket=true,
        test_adjoint_action=true,
        test_diff=true,
    )
end



# TODO: move to GroupTools
function rand_lie(rng::AbstractRNG, G)
    return rand(rng, GroupTools.algebra(G))
end

@testset "eltype rand_lie" begin
    G = MultiAffine(Unitary(4), 3)
    @test eltype(rand_lie(rng, G)) <: Complex
end




"""
``Ad`` and ``exp`` commute:
```math
Ad_{exp(ξ_1)}ξ_2 = exp(ad_{ξ_1}) ξ_2
```
"""
function check_exp_ad(G, vel, tvel)
    χ = exp_lie(G, vel)
    Ad_exp = adjoint_action(G, χ, tvel)

    lie_bracket(G, vel, tvel)
    B = DefaultOrthogonalBasis()
    der = GroupTools.matrix_from_lin_endomorphism(G, ξ -> lie_bracket(G, vel, ξ), B)
    mor = exp(der)
    tvel_coord = get_coordinates_lie(G, tvel, B)
    exp_ad = get_vector_lie(G, mor * tvel_coord, B)
    return isapprox(GroupTools.algebra(G), Ad_exp, exp_ad)
end


@testset "exp (ad_ξ) = Ad_exp(ξ)" for G in [
    MultiDisplacement(3,2),
    MultiDisplacement(2),
    MultiAffine(Unitary(3), 3), # broken: exp(ad_{ξ}) cannot be computed
]
    vel = rand_lie(rng, G)
    tvel = rand_lie(rng, G)
    @test check_exp_ad(G, vel, tvel) broken=G isa MultiAffine{<:Unitary}
end

_adjoint_action(G::MultiAffine, p, X) = begin
    tmp = allocate_result(G, adjoint_action, X)
    return _adjoint_action!(G, tmp, p, X)
end

_adjoint_action!(G::MultiAffine, tmp, p, X) = begin
    mat = affine_matrix(G, p)
    matinv = affine_matrix(G, inv(G, p))
    res = mat * screw_matrix(G, X) * matinv
    map(copyto!, submanifold_components(G, tmp), submanifold_components(G, res))
    return tmp
end

@testset "Compare Adjoint Implementations" for G in [
    MultiDisplacement(3, 2),
    MultiDisplacement(2),
    MultiAffine(Unitary(3), 2),
]
    χ = rand(rng, G)
    ξ = rand_lie(rng, G)
    expected = _adjoint_action(G, χ, ξ)
    computed = adjoint_action(G, χ, ξ)
    @test isapprox(GroupTools.algebra(G), expected, computed)
end

"""
    check_adjoint_action(G, grp_rep, alg_rep, χ, ξ)

The group representation ``ρ`` and algebra representation ``ρ``
 commute with the adjoint action:
```math
ρ(χ) ρ(ξ) ρ(χ)^{-1} = ρ(Ad_χ (ξ))
```
"""
check_adjoint_action(G, grp_rep, alg_rep, χ, ξ) = begin
    mat = grp_rep(G, χ)
    matinv = inv(mat)
    expected = mat * alg_rep(G, ξ) * matinv
    computed = alg_rep(G, adjoint_action(G, χ, ξ))
    return isapprox(computed, expected)
end


"""
    check_inv_rep(G, grp_rep, χ)

The group representation ``ρ`` commutes with the inverse.
```math
ρ(χ)^{-1} = ρ(χ^{-1})
```
"""
check_inv_rep(G, grp_rep, χ) = begin
    computed = grp_rep(G, inv(G, χ))
    expected = inv(grp_rep(G, χ))
    return isapprox(computed, expected)
end



_switch_sign(ξ, ::LeftSide) = ξ
_switch_sign(ξ, ::RightSide) = -ξ
"""
    check_apply_diff_group(G, side::GroupActionSide)

The left group operation action on itself ``α(χ_1)χ_2``
is either (left side)
```math
α(χ_1)χ_2 = χ_1 χ_2
```
or (right side)
```math
α(χ_1)χ_2 = χ_2 χ_1^{-1}
```
Now fix ``χ_2 = 1`` (1 is the identity of ``G``) and define ``f : G → G`` by ``f(χ) := α(χ) 1``. Since ``f(1) = 1``,
its differential at identity is a map ``Alg(G) → Alg(G)``.
This map is either
- ``Id`` (left side)
- ``-Id`` (right side)
"""
check_apply_diff_group(G, ξ, side::Manifolds.GroupActionSide) = begin
    id = identity_element(G)
    ξ_ = apply_diff_group(GroupOperationAction(G, (LeftAction(), side)), id, ξ, id)
    return isapprox(GroupTools.algebra(G), ξ, _switch_sign(ξ_, side))
end

@testset "Diff $G" for G in [
    MultiDisplacement(3, 2),
    SpecialEuclidean(3),
    SpecialOrthogonal(3),
    ]
    @testset "$side" for side in [LeftSide(), RightSide()]
        ξ = rand(rng, GroupTools.algebra(G))
        @test check_apply_diff_group(G, ξ, side)
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
    @testset "MultiAffine(G, $size) creates proper type" for dim in [4]
        for size in [3]
            GA = MultiAffine(Orthogonal(dim), size)
            @test isa(GA, MultiAffine{typeof(Orthogonal(dim)),dim,size,ℝ})
            @test !isa(GA, MultiAffine{typeof(Orthogonal(dim + 1)),dim + 1,size,ℝ})
            @test isa(MultiAffine(Unitary(dim), size), MultiAffine)
        end
    end
end

"""
    check_exp_lie_point(G, ξ)

The Lie group exponential sends the vector ξ
to an element in the group.
"""
check_exp_lie_point(G, ξ) = is_point(G, exp_lie(G, ξ))

"""
    check_adjoint_action_in_alg(G, χ, ξ)

The adjoint action of χ on ξ is an element of Alg(G):
```math
Ad_{χ}ξ ∈ Alg(G)
```
"""
check_adjoint_action_in_alg(G, χ, ξ) = is_vector(G, identity_element(G), adjoint_action(G, χ, ξ))

"""
    check_grp_rep_Identity(G)

The representation works at the Identity(G) point.
"""
check_grp_rep_Identity(G, grp_rep) = begin
    expected = grp_rep(G, identity_element(G))
    computed = grp_rep(G, Identity(G))
    return isapprox(expected, computed)
end

check_zero_Identity(G) = isapprox(GroupTools.algebra(G),
                                  zero_vector(G, Identity(G)),
                                  zero_vector(G, identity_element(G)))

check_from_normal_grp(G::MultiAffine, ts) = begin
    χ1 = from_normal_grp(G, eachcol(ts)...)
    χ2 = from_normal_grp(G, ts)
    return isapprox(G, χ1, χ2)
end

check_from_normal_alg(G::MultiAffine, ts) = begin
    ξ1 = from_normal_alg(G, eachcol(ts)...)
    ξ2 = from_normal_alg(G, ts)
    return isapprox(GroupTools.algebra(G), ξ1, ξ2)
end

"""
    rand_trans(rng, G::MultiAffine)

Random translation part of the group `G`.
"""
rand_trans(rng, G::MultiAffine{TH, dim, size}) where {TH, dim, size} = randn(rng, dim, size)

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
    MultiAffine(Unitary(3), 2),
    ]
    @test check_grp_rep_Identity(G, affine_matrix)
    vel = rand_lie(rng, G)
    pt = rand(rng, G)
    @test check_exp_lie_point(G, vel)
    @test check_adjoint_action_in_alg(G, pt, vel)
    @testset "zero_element" begin
        @test check_zero_Identity(G)
    end
    @testset "Inverse & Matrix" begin
        χ = rand(rng, G)
        @test check_inv_rep(G, affine_matrix, χ)
    end
    @testset "Adjoint action & matrix" begin
        χ = rand(rng, G)
        ξ = rand_lie(rng, G)
        @test check_adjoint_action(G, affine_matrix, screw_matrix, χ, ξ)
    end
    @testset "Lie Bracket & matrix" begin
        v1, v2 = [rand_lie(rng, G) for i in 1:2]
        @test check_alg_rep(G, screw_matrix, v1, v2)
    end
    @testset "Composition & matrix" begin
        p1, p2 = [rand(rng, G) for i in 1:2]
        @test check_grp_rep_compose(G, affine_matrix, p1, p2)
    end
end

"""
    check_grp_rep_compose(G, ρ, χ1, χ2)

The group representation ``ρ`` commutes with composition.
```math
 ρ(χ_1 χ_2) = ρ(χ_1) ρ(χ_2)
```
where the latter is a matrix multiplication.
"""
check_grp_rep_compose(G, grp_rep, χ1, χ2) = begin
    m1, m2 = [grp_rep(G, p) for p in [χ1, χ2]]
    expected = m1 * m2
    computed = grp_rep(G, compose(G, χ1, χ2))
    return isapprox(expected, computed)
end

"""
    check_alg_rep(G, alg_rep, ξ1, ξ2)

The algebra representation ``ρ`` is an algebra morphism.
```math
ρ([ξ_1, ξ_2]) = [ρ(ξ_1), ρ(ξ_2)]
```
where the latter is a matrix commutator.
"""
check_alg_rep(G, alg_rep, ξ1, ξ2) = begin
    m1, m2 = [alg_rep(G, v) for v in [ξ1,ξ2]]
    expected = m1*m2 - m2*m1
    computed = alg_rep(G, lie_bracket(G, ξ1, ξ2))
    return isapprox(expected, computed)
end




include("multiaffine/apply_diff_group.jl")
include("multiaffine/inv_diff.jl")
