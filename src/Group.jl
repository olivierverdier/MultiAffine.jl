



const MultiColumnwiseMultiplicationAction{G,dim,size,𝔽} = Manifolds.ColumnwiseMultiplicationAction{
    LeftAction,
    TranslationGroup{ManifoldsBase.TypeParameter{Tuple{dim,size}},𝔽},
    G,
}

const MultiAffineGroup{G,dim,size,𝔽} = SemidirectProductGroup{
    𝔽,
    TranslationGroup{ManifoldsBase.TypeParameter{Tuple{dim,size}},𝔽},
    G,
    MultiColumnwiseMultiplicationAction{G,dim,size,𝔽}
}

const MultiAffineOp{G,dim,size,𝔽} = Manifolds.SemidirectProductOperation{MultiColumnwiseMultiplicationAction{G,dim,size,𝔽}}

include("Compat.jl")

@doc raw"""
    MultiAffineGroup(G, k=1)

An affine group modelling matrices of the form
```math
χ = \begin{bmatrix}
\mathbf{1} & \mathbf{0} \\
X & g
\end{bmatrix}
```
where ``g`` is a matrix element of the group ``G``, represented in dimension ``n``,
and ``X`` is a ``n × k`` matrix.
If we denote such an element by ``[X,g]``,
the multiplication law is
```math
[X,g] [X',g'] = [X+gX';gg']
```
"""
function MultiAffineGroup(G::Manifolds.GeneralUnitaryMultiplicationGroup{ManifoldsBase.TypeParameter{Tuple{dim}},𝔽}, size::Integer=1) where {dim, 𝔽}
    space = TranslationGroup(dim,size;field=𝔽)
    action = Manifolds.ColumnwiseMultiplicationAction(space, G)
    group = _GroupManifold(ProductManifold(space, G), Manifolds.SemidirectProductOperation(action))
    return group
end


Base.show(io::IO, ::MultiAffineGroup{G, <:Any,size}) where {G, size} = print(io, "MultiAffineGroup($(G), $(size))")

_get_representation_dim(::MultiAffineGroup{<:Any,dim,size}
                        ) where {dim,size} = dim+size


function Manifolds.allocate_result(G::MultiAffineGroup, ::Union{typeof(affine_matrix),typeof(screw_matrix)}, Xis...)
    d = _get_representation_dim(G)
    return allocate(Xis[1], Manifolds.Size(d,d))
end

Base.@propagate_inbounds function Manifolds._padvector!(
    ::MultiAffineGroup{<:Any, <:Any, size},
    X::AbstractMatrix,
) where {size}
    for i in Base.Iterators.take(axes(X, 1), size)
        for j in axes(X, 2)
            X[i, j] = 0
        end
    end
    return X
end

Base.@propagate_inbounds function Manifolds._padpoint!(
    G::MultiAffineGroup{<:Any,<:Any,size},
    q::AbstractMatrix,
) where {size}
    Manifolds._padvector!(G, q)
    for (i,j) in Base.Iterators.take(zip(axes(q)...), size)
        q[i,j] = 1
    end
    return q
end


function Manifolds.affine_matrix(
    G::MultiAffineGroup,
    p
    )
    pis = submanifold_components(G, p)
    pmat = allocate_result(G, affine_matrix, pis...)
    map(copyto!, submanifold_components(G, pmat), pis)
    @inbounds Manifolds._padpoint!(G, pmat)
    return pmat
end

function Manifolds.affine_matrix(
    G::MultiAffineGroup{TH, dim, size, 𝔽},
    ::Identity{MultiAffineOp{TH, dim, size, 𝔽}}
    ) where {TH, dim, size, 𝔽}
    n = _get_representation_dim(G)
    return LinearAlgebra.Diagonal{Float64}(LinearAlgebra.I, n)
end

function Manifolds.zero_vector(
    G::MultiAffineGroup{TH, dim, size, 𝔽},
    ::Identity{MultiAffineOp{TH, dim, size, 𝔽}}
    ) where {TH, dim, size, 𝔽}
    res = allocate_result(G, typeof(zero_vector))
    fill!(res, 0.)
    return res
end


function Manifolds.screw_matrix(G::MultiAffineGroup, X)
    Xis = submanifold_components(G, X)
    Xmat = allocate_result(G, screw_matrix, Xis...)
    map(copyto!, submanifold_components(G, Xmat), Xis)
    @inbounds Manifolds._padvector!(G, Xmat)
    return Xmat
end


Base.@propagate_inbounds function Manifolds.submanifold_component(
    ::MultiAffineGroup{<:Any,dim,size},
    p::AbstractMatrix,
    ::Val{1},
    ) where {dim,size}
    return view(p,
                last(axes(p,1), dim),
                first(axes(p,2), size)
                )
end

Base.@propagate_inbounds function Manifolds.submanifold_component(
    ::MultiAffineGroup{<:Any,dim},
    p::AbstractMatrix,
    ::Val{2},
    ) where {dim}
    return view(p,
                last(axes(p,1), dim),
                last(axes(p,2), dim)
                )
end

function Manifolds.submanifold_components(
    G::MultiAffineGroup,
    p::AbstractMatrix,
    )
    d = _get_representation_dim(G)
    @assert Base.size(p) == (d,d)
    @inbounds t = submanifold_component(G, p, Val(1))
    @inbounds R = submanifold_component(G, p, Val(2))
    return (t, R)
end


function Manifolds._log_lie!(G::MultiAffineGroup, X, q)
    qmat = affine_matrix(G, q)
    Xmat = Manifolds.log_safe(qmat)
    map(copyto!, submanifold_components(G, X), submanifold_components(G, Xmat))
    Manifolds._padvector!(G, X)
    return X
end

function Manifolds.exp_lie!(G::MultiAffineGroup, q, X)
    Xmat = screw_matrix(G, X)
    qmat = exp(Xmat)
    map(copyto!, submanifold_components(G, q), submanifold_components(G, qmat))
    Manifolds._padpoint!(G, q)
    return q
end


Base.exp(G::MultiAffineGroup, χ, v) = ManifoldGroupUtils.exp_group(G, χ, v)
Manifolds.exp!(G::MultiAffineGroup, tmp, χ, v) = ManifoldGroupUtils.exp_group!(G, tmp, χ, v)
Base.log(G::MultiAffineGroup, χ1, χ2) = ManifoldGroupUtils.log_group(G, χ1, χ2)
Manifolds.log!(G::MultiAffineGroup, tmp, χ1, χ2) = ManifoldGroupUtils.log_group!(G, tmp, χ1, χ2)


Manifolds.adjoint_action!(G::MultiAffineGroup, tmp, p, X, ::LeftAction) = begin
    np, hp = submanifold_components(G, p)
    n, h = submanifold_components(G, tmp)
    nX, hX = submanifold_components(G, X)
    H = factor_group(G)
    adjoint_action!(H, h, hp, hX)
    A = G.op.action
    apply!(A, n, hp, nX)
    LinearAlgebra.axpy!(-1, apply_diff_group(A, Identity(H), h, np), n)
    return tmp
end





function Manifolds.lie_bracket(G::MultiAffineGroup, v1, v2, )
    res = allocate_result(G, lie_bracket)
    return lie_bracket!(G, res, v1, v2)
end

function Manifolds.lie_bracket!(G::MultiAffineGroup, res, v1, v2, )
    A = G.op.action
    n1, h1 = submanifold_components(v1)
    n2, h2 = submanifold_components(v2)
    rn = apply(A, h1, n2) - apply(A, h2, n1)
    rh = lie_bracket(factor_group(G), h1, h2)
    map(copyto!, submanifold_components(G, res), [rn,rh])
    return res
end

# The two morphisms corresponding to the sequence
# 0 -> V -> G -> H -> 0

function _fill_in!(G::MultiAffineGroup, x, ts)
    X = submanifold_component(G, x, 1)
    copyto!(X, ts)
    return x
end

_fill_in!(G::MultiAffineGroup, x, ts...) = _fill_in!(G, x, hcat(ts...))

factor_group(G::MultiAffineGroup) = submanifold(G, 2)
normal_group(G::MultiAffineGroup) = submanifold(G, 1)

from_normal_grp(G::MultiAffineGroup, ts...) = begin
    x = identity_element(G)
    return _fill_in!(G, x, ts...)
end

from_normal_alg(G::MultiAffineGroup, ts...) = begin
    x = zero_vector(G, identity_element(G))
    return _fill_in!(G, x, ts...)
end

to_factor(G::MultiAffineGroup, pt) = submanifold_component(G, pt, 2)
to_factor_grp(G::MultiAffineGroup, χ) = to_factor(G, χ)
to_factor_alg(G::MultiAffineGroup, ξ) = to_factor(G, ξ)

to_normal(G::MultiAffineGroup, pt) = submanifold_component(G, pt, 1)
to_normal_grp(G::MultiAffineGroup, χ) = to_normal(G, χ)
to_normal_alg(G::MultiAffineGroup, ξ) = to_normal(G, ξ)

"""
    normal_indices(::MultiAffineGroup, idx)

If ``[X,R]`` is an element of the group ``G``,
return the indices of the normal component `X`.
"""
normal_indices(::MultiAffineGroup{<:Any,dim,size}, idx) where {dim,size} = collect(Iterators.take(idx, dim * size))

"""
    normal_indices_at(G::MultiAffine, idx, pos=0)

If ``[X,R]`` is an element of the group ``G``,
return the indices of the column vector of ``X`` at offset `pos` (starting at zero).
"""
normal_indices_at(G::MultiAffineGroup{<:Any,dim}, idx, pos=0) where {dim} = collect(Iterators.take(Iterators.drop(normal_indices(G, idx), pos * dim), dim))

"""
    factor_indices(G::MultiAffineGroup{<:Any, dim, size}, idx)

If ``[X,R]`` is an element of the group ``G``,
return the indices of the factor component `R`.
"""
factor_indices(::MultiAffineGroup{<:Any, dim, size}, idx) where {dim,size} = collect(Iterators.drop(idx, size * dim))
