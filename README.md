# MultiAffine

[![Build Status](https://github.com/olivierverdier/MultiAffine.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/olivierverdier/MultiAffine.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/olivierverdier/MultiAffine.jl/graph/badge.svg?token=aTe2GSxvIw)](https://codecov.io/gh/olivierverdier/MultiAffine.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierverdier.github.io/MultiAffine.jl/)


Given a group $`H`$, represented in dimension $`n`$ over some field $`\mathbf{F}`$
this package implements the group
```math
(\mathbf{F}^n)^k \rtimes H
```
The group consists of matrices of the form
```math
[X;h] \equiv \begin{bmatrix}
\mathbf{1} & \mathbf{0} \\
X & h
\end{bmatrix}
```
where $`h`$ is a matrix belonging to the group $`H`$,
and $`X`$ is a $`n × k`$ matrix.
If we denote such an element by $`[X;h]`$,
the multiplication law is
```math
[X;h] [X';h'] = [X+hX';hh']
```

For instance, to implement the group
```math
(\mathbf{C}^2)^2 \rtimes U(2)
```
you can use
```julia
using MultiAffine
using Manifolds
G = MultiAffineGroup(Unitary(2), 2)
identity_element(G) # ([0.0 0.0; 0.0 0.0], ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im])
```

When the group $`H`$ is the [special orthogonal group](https://en.wikipedia.org/wiki/Orthogonal_group) $`SO(n)`$, one can use the alias `MultiDisplacementGroup(n,k)` to implement the group
```math
(\mathbf{R}^n)^k \rtimes SO(n)
```
This has the same effect as calling `MultiAffineGroup(SpecialOrthogonal(n), k)`. 

```julia
G = MultiDisplacementGroup(3,2)
identity_element(G) # ([0.0 0.0; 0.0 0.0; 0.0 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
```
