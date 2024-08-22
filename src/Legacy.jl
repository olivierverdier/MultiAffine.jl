"""
This file contains code which will appear in Manifolds.jl v0.10
and is only needed for this package to work with older versions of Manifolds.jl.
"""
function Manifolds.adjoint_action(G::MultiAffineGroup, p, X, conv=LeftAction())
    tmp = allocate_result(G, adjoint_action, X)
    return adjoint_action!(G, tmp, p, X, conv)
end

Manifolds.adjoint_action(::MultiAffineGroup, ::Identity, X, ::LeftAction) = X
Manifolds.adjoint_action(::MultiAffineGroup, ::Identity, X, ::RightAction) = X
Manifolds.adjoint_action!(G::MultiAffineGroup, Y, ::Identity, X, ::LeftAction) = copyto!(G, Y, X)

inverse_adjoint_action!(G::MultiAffineGroup, Y, p, X, conv=LeftAction()) = adjoint_action!(G, Y, inv(G, p), X, conv)

Manifolds.translate_diff!(G::MultiAffineGroup,
    Y, ::Any, ::Any, X,
    ::Manifolds.LeftForwardAction,
) = copyto!(G, Y, X)
Manifolds.translate_diff!(G::MultiAffineGroup,
    Y, ::Any, ::Any, X,
    ::Manifolds.RightForwardAction,
) = copyto!(G, Y, X)
Manifolds.translate_diff!(
    G::MultiAffineGroup,
    Y, p, ::Any, X,
    ::Manifolds.LeftBackwardAction,
) = adjoint_action!(G, Y, p, X)
Manifolds.translate_diff!(
    G::MultiAffineGroup,
    Y, p, ::Any, X,
    ::Manifolds.RightBackwardAction,
) = inverse_adjoint_action!(G, Y, p, X)


Manifolds.inv_diff(::MultiAffineGroup, ::Identity, ξ) = -ξ
Manifolds.inv_diff!(G::MultiAffineGroup, Y, p, X) = begin
    adjoint_action!(G, Y, p, X)
    Y .*= -1
    return Y
end


