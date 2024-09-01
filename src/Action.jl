"""
    MultiAffineAction(
      group::MultiAffineGroup,
      selector::AbstractArray,
      conv::ActionDirection=LeftAction()
      )

Given a fixed vector ``S`` of size ``k`` (the `selector`) (or a matrix of size ``k×m``),
this defines an action of the element ``[X;R]`` of the [`MultiAffineGroup`](@ref) group (so ``X`` is a ``n×k`` matrix and ``R`` is an element of a matrix group)
 on the vector ``p`` of size ``n`` (or the matrix ``p`` of size ``n×m``).
The action is defined by
```math
[X;R]⋅p := XS+Rp
```
"""
struct MultiAffineAction{TAD<:ActionDirection,TG,TS<:AbstractArray} <: AbstractGroupAction{TAD}
    group::TG
    selector::TS # vector of length `size`
end

Base.show(io::IO, A::MultiAffineAction{TAD}) where {TAD} = print(io, "MultiAffineAction($(A.group), $(A.selector), $TAD())")

function MultiAffineAction(
    group::MultiAffineGroup{<:Any, <:Any, size},
    selector::AbstractArray,
    conv::ActionDirection=LeftAction()
    ) where {size}
    @assert Base.size(selector, 1) == size
    return MultiAffineAction{typeof(conv), typeof(group), typeof(selector)}(group, selector)
end

"""
    MultiAffineAction(
      group::MultiAffineGroup{<:Any, <:Any, 1},
      conv::ActionDirection=LeftAction()
      )

The standard affine action, obtained with a selector equal to ``[1]``.
"""
function MultiAffineAction(
    group::MultiAffineGroup{<:Any, <:Any, 1},
    conv::ActionDirection=LeftAction()
    )
    return MultiAffineAction(group, [1], conv)
end

function Manifolds.switch_direction(A::MultiAffineAction{TAD}) where {TAD}
    return MultiAffineAction(A.group, A.selector, switch_direction(TAD()))
end


Manifolds.base_group(A::MultiAffineAction) = A.group
Manifolds.group_manifold(::MultiAffineAction{<:Any,<:MultiAffineGroup{<:Any, dim, <:Any, 𝔽}, <:AbstractVector}) where {dim,𝔽} = Euclidean(dim; field=𝔽)
Manifolds.group_manifold(A::MultiAffineAction{<:Any,<:MultiAffineGroup{<:Any, dim, <:Any, 𝔽}, <:AbstractMatrix}) where {dim,𝔽} = Euclidean(dim, size(A.selector, 2); field=𝔽)

get_selector(A::MultiAffineAction) = A.selector

Manifolds.apply(::MultiAffineAction{<:Any, <:MultiAffineGroup{TH,dim,size,𝔽}},
                ::Identity{MultiAffineOp{TH,dim,size,𝔽}}, p) where {TH,dim,size,𝔽} = p

function Manifolds.apply!(
    A::MultiAffineAction{LeftAction},
    q,
    χ,
    p, # vector of size dim
    )
    G = base_group(A)
    M,R = submanifold_components(G, χ)
    sel = get_selector(A)
    # compute: M*sel + R*p
    LinearAlgebra.mul!(q, M, sel)
    LinearAlgebra.mul!(q, R, p, 1, 1)
    return q
end


Manifolds.apply!(A::MultiAffineAction{RightAction}, q, χ, p) = apply!(switch_direction(A), q, inv(base_group(A), χ), p)



function Manifolds.apply_diff_group(
    A::MultiAffineAction{LeftAction, <:MultiAffineGroup{TH,dim,size,𝔽}},
    ::Identity{MultiAffineOp{TH,dim,size,𝔽}},
    ξ,
    p
    ) where {TH,dim,size,𝔽}
    G = base_group(A)
    μ, ω = submanifold_components(G, ξ)
    sel = get_selector(A)
    @assert is_point(group_manifold(A), p)
    return μ*sel + ω*p
end


function Manifolds.apply_diff_group(
    A::MultiAffineAction{RightAction, <:MultiAffineGroup{TH,dim,size,𝔽}},
    I::Identity{MultiAffineOp{TH,dim,size,𝔽}},
    ξ,
    p
    ) where {TH,dim,size,𝔽}
    return apply_diff_group(switch_direction(A), I, -ξ, p)
end
