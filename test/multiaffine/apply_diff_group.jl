using Test
using Manifolds
import GroupTools
import Random
rng = Random.default_rng()



test_apply_diff_group_set(G, χ, ξ, χ_) = begin
    @testset "$dir, $side" for side in [LeftSide(), RightSide()], dir in [LeftAction(), RightAction()]
            @test GroupTools.check_apply_diff_group(GroupOperationAction(G, (dir, side)), χ, ξ, χ_)
    end
end


@testset "apply_diff_group" begin
    G = MultiDisplacement(3,2)
    # χ1 = rand(rng, G)
    χ1 = Identity(G)
    χ2 = rand(rng, G)
    ξ = rand(rng, TangentSpace(G, identity_element(G)))
    test_apply_diff_group_set(G, χ1, ξ, χ2)
end
