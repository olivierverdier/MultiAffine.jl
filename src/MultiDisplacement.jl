
const MultiDisplacementGroup{dim,size} = MultiAffineGroup{
    SpecialOrthogonal{ManifoldsBase.TypeParameter{Tuple{dim}}},
    dim,
    size,
    ℝ
}

@deprecate MultiDisplacement(dim::Integer, size::Integer=1) MultiDisplacementGroup(dim, size)
@deprecate MultiDisplacement  MultiDisplacementGroup

"""
    MultiDisplacementGroup(n, k=1)

A special case of the [`MultiAffineGroup`](@ref) group, where the
underlying group is the special orthogonal group ``SO(n)``.
This is just a convenient alias
```julia
MultiDisplacementGroup(n, k=1) = MultiAffineGroup(SpecialOrthogonal(n), k)
```

When ``k=1`` (the default), this is the [*Special Euclidean Group*](https://en.wikipedia.org/wiki/Euclidean_group).
"""
MultiDisplacementGroup(dim::Integer, size::Integer=1) = MultiAffineGroup(SpecialOrthogonal(dim), size)

Base.show(io::IO, ::MultiDisplacementGroup{dim,size}) where {dim,size} = print(io, "MultiDisplacementGroup($(dim), $(size))")

Base.show(io::IO, ::MultiDisplacementGroup{dim,1}) where {dim} = print(io, "MultiDisplacementGroup($(dim))")
