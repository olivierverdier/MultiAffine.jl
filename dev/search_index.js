var documenterSearchIndex = {"docs":
[{"location":"#MultiAffine.jl","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"","category":"section"},{"location":"","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"This package consists of an implementation of","category":"page"},{"location":"","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"the MultiAffine group\nan action on vector spaces","category":"page"},{"location":"#MultiAffine-Group","page":"MultiAffine.jl","title":"MultiAffine Group","text":"","category":"section"},{"location":"","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"Modules = [MultiAffine]\nOrder = [:type, :function]\nPages = [\"Group.jl\"]","category":"page"},{"location":"#MultiAffine.MultiAffineGroup-Union{Tuple{Manifolds.GeneralUnitaryMultiplicationGroup{ManifoldsBase.TypeParameter{Tuple{dim}}, 𝔽, S} where S}, Tuple{𝔽}, Tuple{dim}, Tuple{Manifolds.GeneralUnitaryMultiplicationGroup{ManifoldsBase.TypeParameter{Tuple{dim}}, 𝔽, S} where S, Integer}} where {dim, 𝔽}","page":"MultiAffine.jl","title":"MultiAffine.MultiAffineGroup","text":"MultiAffineGroup(G, k=1)\n\nAn affine group modelling matrices of the form\n\nχ = beginbmatrix\nmathbf1  mathbf0 \nX  g\nendbmatrix\n\nwhere g is a matrix element of the group G, represented in dimension n, and X is a n  k matrix. If we denote such an element by Xg, the multiplication law is\n\nXg Xg = X+gXgg\n\n\n\n\n\n","category":"method"},{"location":"","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"Modules = [MultiAffine]\nOrder = [:type, :function]\nPages = [\"MultiDisplacement.jl\"]","category":"page"},{"location":"#MultiAffine.MultiDisplacement","page":"MultiAffine.jl","title":"MultiAffine.MultiDisplacement","text":"MultiDisplacement(n, k=1)\n\nA special case of the MultiAffineGroup group, where the underlying group is the special orthogonal group SO(n). This is just a convenient alias\n\nMultiDisplacement(n, k=1) = MultiAffineGroup(SpecialOrthogonal(n), k)\n\nWhen k=1 (the default), this is the Special Euclidean Group.\n\n\n\n\n\n","category":"type"},{"location":"#MultiAffine-Action","page":"MultiAffine.jl","title":"MultiAffine Action","text":"","category":"section"},{"location":"","page":"MultiAffine.jl","title":"MultiAffine.jl","text":"Modules = [MultiAffine]\nOrder = [:type, :function]\nPages = [\"Action.jl\"]","category":"page"},{"location":"#MultiAffine.MultiAffineAction","page":"MultiAffine.jl","title":"MultiAffine.MultiAffineAction","text":"MultiAffineAction(\n  group::MultiAffineGroup,\n  selector::AbstractArray,\n  conv::ActionDirection=LeftAction()\n  )\n\nGiven a fixed vector S of size k (the selector) (or a matrix of size km), this defines an action of the element XR of the MultiAffineGroup group (so X is a nk matrix and R is an element of a matrix group)  on the vector p of size n (or the matrix p of size nm). The action is defined by\n\nXRp = XS+Rp\n\n\n\n\n\n","category":"type"},{"location":"#MultiAffine.MultiAffineAction-2","page":"MultiAffine.jl","title":"MultiAffine.MultiAffineAction","text":"MultiAffineAction(\n  group::MultiAffineGroup{<:Any, <:Any, 1},\n  conv::ActionDirection=LeftAction()\n  )\n\nThe standard affine action, obtained with a selector equal to 1.\n\n\n\n\n\n","category":"type"}]
}
