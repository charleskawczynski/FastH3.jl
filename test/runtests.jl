using Test

include(joinpath(@__DIR__, "..", "src", "H3X.jl"))
using .H3X

@testset "H3X" begin
    include("test_layer0.jl")
    include("test_layer1_2.jl")
    include("test_layer3_4.jl")
end
