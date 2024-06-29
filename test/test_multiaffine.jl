using Test
using MultiAffinity

import GroupTools

import ManifoldsBase as MB
using Manifolds
import Random: rand!, AbstractRNG

import Random

rng = Random.default_rng()



# TODO: this is also implemented in Motion.jl
get_adjoint_matrix(G, vel, B::AbstractBasis) = GroupTools.matrix_from_lin_endomorphism(G, ξ -> lie_bracket(G, vel, ξ), B)


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
    end
    @testset "MultiAffine(G, size) creates proper type" begin
        GA = MultiAffine(Orthogonal(4), 3)
        @test isa(GA, MultiAffine{typeof(Orthogonal(4)), 4, 3, ℝ})
        @test !isa(GA, MultiAffine{typeof(Orthogonal(5)), 5, 3, ℝ})
        @test isa(MultiAffine(Unitary(4), 5), MultiAffine)
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

end

function test_multi_affine(rng, G::MultiAffine{TH,dim,size,𝔽}
                           ) where {TH,dim,size,𝔽}
  @testset "Test $(repr(G))" begin
      vel = rand_lie(rng, G)
      pt = rand(rng, G)
      x = exp_lie(G, vel)
      @test is_point(G, x)
      v_ = adjoint_action(G, pt, vel)
      @test is_vector(G, identity_element(G), v_)
      affine_matrix(G, Identity(G))
      @testset "zero_element" begin
          z = zero_vector(G, Identity(G))
          z_ = zero_vector(G, identity_element(G))
          @test isapprox(G, z, z_)
      end
      @testset "from/to" begin
          ts = randn(rng, dim, size)
          χ1 = from_normal_grp(G, eachcol(ts)...)
          χ2 = from_normal_grp(G, ts)
          @test isapprox(G, χ1, χ2)
          ξ1 = from_normal_alg(G, eachcol(ts)...)
          ξ2 = from_normal_alg(G, ts)
          @test isapprox(G, Identity(G), ξ1, ξ2)
      end
      @testset "Lie Bracket & matrix" begin
          v1, v2 = [rand_lie(rng, G) for i in 1:2]
          m1, m2 = [screw_matrix(G, v) for v in [v1,v2]]
          comm = m1*m2 - m2*m1
          expected = ArrayPartition(submanifold_components(G, comm)...)
          computed = lie_bracket(G, v1, v2)
          @test isapprox(G, Identity(G),  expected, computed)
      end
      @testset "Composition & matrix" begin
          p1, p2 = [rand(rng, G) for i in 1:2]
          m1, m2 = [affine_matrix(G, p) for p in [p1,p2]]
          prod = m1*m2
          expected = ArrayPartition(submanifold_components(G, prod)...)
          computed = compose(G, p1, p2)
          @test isapprox(G, expected, computed)
      end
  end
end



test_exp_ad(Random.default_rng(), MultiDisplacement(3, 2))

@testset "MultiAffine" for G in
    [
        MultiDisplacement(3, 2),
        # MultiAffine(Unitary(4), 3),
    ]
    # begin
    test_multi_affine(Random.default_rng(), G)
    @testset "Adjoint action & matrix" begin
        χ = rand(rng, G)
        ξ = rand_lie(rng, G)
        @test check_adjoint_action(G, affine_matrix, screw_matrix, χ, ξ)
    end
end

include("multiaffine/apply_diff_group.jl")
include("multiaffine/inv_diff.jl")
