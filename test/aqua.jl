import Aqua
import SamApp2025
using Test: @testset

@testset verbose = true "aqua" begin
    # persistent_tasks is disabled because SamApp2025 boots a Python interpreter
    # at load time (Logomaker -> PythonCall), which keeps Aqua's helper
    # precompilation process from exiting cleanly.
    Aqua.test_all(SamApp2025; ambiguities = false, persistent_tasks = false)
end
