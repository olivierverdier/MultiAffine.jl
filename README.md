# MultiAffine

[![Build Status](https://github.com/olivierverdier/MultiAffine.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/olivierverdier/MultiAffine.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/olivierverdier/MultiAffine.jl/graph/badge.svg?token=aTe2GSxvIw)](https://codecov.io/gh/olivierverdier/MultiAffine.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierverdier.github.io/MultiAffine.jl/)


This package models an affine group consisting of matrices of the form
```math
χ = \begin{bmatrix}
\mathbf{1} & \mathbf{0} \\
X & h
\end{bmatrix}
```
where $`h`$ is a matrix element of some group $`H`$, represented in dimension $`n`$,
and $`X`$ is a $`n × k`$ matrix.
If we denote such an element by $`[X,h]`$,
the multiplication law is
```math
[X,h] [X',h'] = [X+hX';hh']
```

One example of such group is
```julia
using MultiAffine
using Manifolds
G = MultiAffineGroup(Unitary(2), 2)
identity_element(G) # ([0.0 0.0; 0.0 0.0], ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im])
```

When the group $`H`$ is the [special orthogonal group](https://en.wikipedia.org/wiki/Orthogonal_group), one can use the alias `MultiDisplacement(n,k)` to create the group `MultiAffine(SpecialOrthogonal(n), k)`, for instance:
```julia
G = MultiDisplacement(3,2)
identity_element(G) # ([0.0 0.0; 0.0 0.0; 0.0 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
```
