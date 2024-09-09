using Test
using Manifolds
import Random
rng = Random.default_rng()



@testset "apply_diff_group $G" for
    G in [MultiDisplacementGroup(3, 2), MultiDisplacementGroup(2)]
    χ1 = rand(rng, G)
    χ2 = rand(rng, G)
    ξ = rand(rng, TangentSpace(G, identity_element(G)))
    @testset "$χ $dir $side" for
        χ in [χ1, Identity(G)],
        side in [LeftSide(), RightSide()],
        dir in [LeftAction(), RightAction()]
            @test GT.check_apply_diff_group(GroupOperationAction(G, (dir, side)), χ1, ξ, χ2)
    end
end
