import Aqua
import SamApp2025
using Test: @testset

@testset verbose = true "aqua" begin
    Aqua.test_all(SamApp2025; ambiguities = false)
end
