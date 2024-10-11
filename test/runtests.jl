include("header.jl")

using MultiAffine

import ManifoldGroupTesting as GT

test_files = [
    "test_group.jl",
    "test_action.jl",
]

include_tests(test_files)
