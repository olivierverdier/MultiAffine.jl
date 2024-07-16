module MultiAffine

# MultiAffine
export MultiAffineGroup,
    MultiDisplacement,
    from_normal_grp, from_normal_alg,
    to_factor_grp, to_factor_alg,
    to_normal_grp, to_normal_alg,
    normal_group, factor_group

export MultiAffineAction, get_selector

import Manifolds
import ManifoldsBase
import Manifolds: # Manifolds
    ℝ,
    ProductManifold, submanifold_component, submanifold_components, submanifold,
    allocate_result, allocate,
    zero_vector,
    is_point
import Manifolds: # Group
    GroupManifold,
    Identity, identity_element,
    TranslationGroup, Euclidean,
    SemidirectProductGroup, affine_matrix, screw_matrix,
    SpecialOrthogonal,
    lie_bracket, lie_bracket!
import Manifolds: # group actions
    apply, apply!, apply_diff_group,
    base_group, group_manifold,
    AbstractGroupAction,
    ActionDirection, LeftAction, RightAction, switch_direction, #--
    adjoint_action, adjoint_action!
import Manifolds:
    RotationActionOnVector, RotationAction

import LinearAlgebra


include("Group.jl")
include("MultiDisplacement.jl")

include("Action.jl")


end
