module MultiAffine

# MultiAffine
export MultiAffineGroup,
    MultiDisplacement,
    from_normal_grp, from_normal_alg,
    to_factor_grp, to_factor_alg,
    to_normal_grp, to_normal_alg,
    normal_group, factor_group

export MultiAffineAction, get_selector
import RecursiveArrayTools

import ManifoldsBase
using Manifolds

import LinearAlgebra


include("Group.jl")
include("MultiDisplacement.jl")

include("Action.jl")


end
