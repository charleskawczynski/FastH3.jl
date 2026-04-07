using Test
using FastH3

@testset "FastH3" begin
    include("test_layer0.jl")
    include("test_layer1_2.jl")
    include("test_layer3_4.jl")
    include("test_callbacks.jl")
    include("test_quality.jl")
end
