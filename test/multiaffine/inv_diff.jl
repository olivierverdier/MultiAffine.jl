using Manifolds
import Random
rng = Random.default_rng()





@testset "inv_diff" begin
    G = MultiDisplacementGroup(3, 2)
    χ = rand(rng, G)
    ξ = rand(rng, TangentSpace(G, identity_element(G)))
    @testset "$side" for side in [LeftSide(), RightSide()]
        @test GT.check_inv_diff(G, χ, ξ, side)
    end
end
