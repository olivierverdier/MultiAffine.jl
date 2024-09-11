using MultiAffine

using Manifolds
import Random: default_rng

rng = default_rng()


"""
    random_multiaffine_action(rng, G::MultiAffineGroup, conv)

Given a multiaffine group,
create a multiaffine action defined by a random selector.
"""
random_multiaffine_action(rng, G::MultiAffineGroup{<:Any, <:Any, size}, conv) where {size} = begin
    sel = randn(rng, size)
    return MultiAffineAction(G, sel, conv)
end


@testset "Test action morphism $G $conv" for G in [
    MultiDisplacementGroup(4,2),
    MultiDisplacementGroup(4),
    ], conv in [LeftAction(), RightAction()]
    A = random_multiaffine_action(rng, G, conv)
    χ1 = rand(rng, G)
    χ2 = rand(rng, G)
    p = rand(rng, group_manifold(A))
    @test GT.check_action_morphism(A, χ1, χ2, p)
end


check_repr(A::MultiAffineAction{conv}) where {conv} = begin
    G = base_group(A)
    sel = get_selector(A)
    expected = "MultiAffineAction($G, $sel, $conv())"
    computed = repr(A)
    return computed == expected
end

"""
    check_default_selector(G::MultiAffineGroup)

If `size` is one, the group G is some
displacement group, and the multiaffine action
can be created with `MultiAffineAction(G)`
with a default selector equal to `[1]`.
"""
check_default_selector(G::MultiAffineGroup{<:Any, <:Any, 1}) = begin
    A = MultiAffineAction(G)
    return get_selector(A) == [1]
end


get_size(::MultiAffineGroup{<:Any,<:Any,size}) where {size} = size
get_dim(::MultiAffineGroup{<:Any,dim}) where {dim} = dim


inverse_if(::Any, χ, ::LeftAction) = χ
inverse_if(G, χ, ::RightAction) = inv(G, χ)

"""
    check_selector(G::MultiDisplacementGroup, χ, k, conv)

For a multiaffine action ``α`` with selector
```math
s = [0,...,1,...0]
```
with a one at the ``k``th place, one should have
```math
α([M,R], 0) = M[:,k]
```
that is, the action on zero “selects”
the ``k``th column of ``M``.
"""
check_selector(G::MultiDisplacementGroup{<:Any,size}, χ, k, conv) where {size} = begin
    sel = zeros(size)
    sel[k] = 1.0
    A = MultiAffineAction(G, sel, conv)
    dim = manifold_dimension(group_manifold(A))
    p = zeros(dim)
    computed = apply(A, χ, p)
    expected = to_normal_grp(G, inverse_if(G, χ, conv))[:, k]
    return isapprox(group_manifold(A), computed, expected)
end


@testset "Test Multiaffine Action" for G in [MultiDisplacementGroup(3,2)]
    for conv in [LeftAction(), RightAction()]
        @testset "Simple Selector" begin
            k = rand(rng, eachindex(zeros(get_size(G))))
            @test check_selector(G, rand(rng, G), k, conv)
        end
        @testset "repr" begin
            A = random_multiaffine_action(rng, G, conv)
            @test check_repr(A)
        end

        @testset "Apply identity" begin
            A = random_multiaffine_action(rng, G, conv)
            p = rand(rng, group_manifold(A))
            GT.check_apply_morphism_Identity(A, p)
        end
        @testset "Trivial action" begin
            A = random_multiaffine_action(rng, G, conv)
            p = rand(rng, group_manifold(A))
            @test GT.check_trivial_infinitesimal_action(A, p)
        end
    end
end

@testset "default selector" for dim in [3]
    se = MultiDisplacementGroup(dim, 1)
    @test check_default_selector(se)
end

"""
    check_multiaffine_action(α::MultiAffineAction, χ, p)

The multiafine action ``α`` with selector ``σ`` satisfies
for any ``χ = [M; R]`` and point ``p``:
```math
α([M;R], p) = Mσ + Rp
```
"""
check_multiaffine_action(A::MultiAffineAction, χ, p) = begin
    computed = apply(A, χ, p)
    M, R = submanifold_components(base_group(A), χ)
    sel = get_selector(A)
    expected = M * sel + R * p
    return isapprox(group_manifold(A), expected, computed)
end


@testset "MultiAffineAction apply $G×$prod" for G in [
    MultiDisplacementGroup(3, 2),
    MultiDisplacementGroup(2)
    ]
    for prod in [1, 2]
        size = get_size(G)
        dim = get_dim(G)
        @testset "Product $sel" for sel in [randn(rng, size), randn(rng, size, prod)]
            A = MultiAffineAction(G, sel)
            @testset "Correct manifold" begin
                expected_manifold = ndims(sel) == 1 ? Euclidean(dim) : Euclidean(dim, prod)
                @test group_manifold(A) == expected_manifold
            end
            @testset "Action" begin
                χ = rand(rng, G)
                p = rand(rng, group_manifold(A))
                @test check_multiaffine_action(A, χ, p)
            end
            @testset "Switch direction $conv" for conv in [LeftAction(), RightAction()]
                A = MultiAffineAction(G, sel, conv)
                χ = rand(rng, G)
                p = rand(rng, group_manifold(A))
                @test GT.check_switch_action_direction(A, χ, p)
            end
        end
    end
end

@testset "MultiAffineAction(group, selector)" begin
    G = MultiDisplacementGroup(3, 2)
    if VERSION >= v"1.8"
        @test_throws ["MethodError", "MultiAffineAction"] MultiAffineAction(G, RightAction())
    end
end
