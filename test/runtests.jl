using MultiAffinity
using Test

test_files = [
    "test_multiaffine.jl",
    "test_multiaffineaction.jl",
]

@time @testset " " for path in test_files
    printstyled("─"^16 * "[ $path ]\n"; color=:yellow)
    @time include(path)
end

