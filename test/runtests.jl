include("header.jl")

using MultiAffine

import ManifoldGroupTesting as GT

test_files = [
    "test_group.jl",
    "test_action.jl",
]

@time @testset " " for path in test_files
    printstyled("─"^16 * "[ $path ]\n"; color=:yellow)
    @time include(path)
end

