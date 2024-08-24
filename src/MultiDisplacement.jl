
const MultiDisplacement{dim,size} = MultiAffineGroup{
    SpecialOrthogonal{ManifoldsBase.TypeParameter{Tuple{dim}}},
    dim,
    size,
    ℝ
}

"""
    MultiDisplacement(n, k=1)

A special case of the [`MultiAffineGroup`](@ref) group, where the
underlying group is the special orthogonal group ``SO(n)``.
This is just a convenient alias
```julia
MultiDisplacement(n, k=1) = MultiAffineGroup(SpecialOrthogonal(n), k)
```

When ``k=1`` (the default), this is the [*Special Euclidean Group*](https://en.wikipedia.org/wiki/Euclidean_group).
"""
MultiDisplacement(dim::Integer, size::Integer=1) = MultiAffineGroup(SpecialOrthogonal(dim), size)

Base.show(io::IO, ::MultiDisplacement{dim,size}) where {dim,size} = print(io, "MultiDisplacement($(dim), $(size))")

Base.show(io::IO, ::MultiDisplacement{dim,1}) where {dim} = print(io, "MultiDisplacement($(dim))")
