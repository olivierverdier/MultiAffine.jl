using Test
using Manifolds
import Random
import GroupTools
rng = Random.default_rng()



test_inv_diff_set(G, χ, ξ) = @testset "inv_diff" for side in [LeftSide(), RightSide()]
    @test GroupTools.check_inv_diff(G, χ, ξ, side)
end


@testset "inv_diff" begin
    G = MultiDisplacement(3,2)
    χ = rand(rng, G)
    ξ = rand(rng, TangentSpace(G, identity_element(G)))
    test_inv_diff_set(G, χ, ξ)
end
