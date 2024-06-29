using Test
using Manifolds
import GroupTools
import Random
rng = Random.default_rng()



@testset "apply_diff_group" begin
    G = MultiDisplacement(3, 2)
    # χ1 = rand(rng, G)
    χ1 = Identity(G)
    χ2 = rand(rng, G)
    ξ = rand(rng, TangentSpace(G, identity_element(G)))
    @testset "$dir, $side" for side in [LeftSide(), RightSide()], dir in [LeftAction(), RightAction()]
            @test GroupTools.check_apply_diff_group(GroupOperationAction(G, (dir, side)), χ1, ξ, χ2)
    end
end
