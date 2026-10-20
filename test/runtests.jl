using DrWatson, Test
# This repo is a DrWatson script project (no src/FastDysonSolver.jl module):
# activate the parent project explicitly, since @quickactivate from test/
# would find the nameless test project instead.
using Pkg
Pkg.activate(dirname(@__DIR__))

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
include(joinpath(@__DIR__, "..", "src", "Junction.jl"))

println("Starting tests")
ti = time()

@testset "FastDysonSolver tests" begin
    @test 1 == 1

    @testset "smoke: junction run on small grid" begin
        p = JunctionParameters(δt = 0.5, T = 5.0, Γl = 1.0, Γr = 1.0,
                               β = 10.0, Δ = 1.0)
        result = simulate_junction(p)
        @test size(matrix(result[:G_R])) == (2 * N(p), 2 * N(p))
        current = compute_average_current(result)
        @test length(current) == N(p)
        @test all(isfinite, current)
    end

    @testset "compat: backend interface (PR 54)" begin
        if isdefined(NonEquilibriumGreenFunction, :with_thread_budget)
            @test NonEquilibriumGreenFunction.current_thread_budget() > 0
            I_budgeted = with_thread_budget(Threads.nthreads()) do
                p = JunctionParameters(δt = 0.5, T = 5.0, Γl = 1.0, Γr = 1.0,
                                       β = 10.0, Δ = 1.0)
                single_run_current(p).I
            end
            @test all(isfinite, I_budgeted)
        else
            @info "NonEquilibriumGreenFunction 0.2: skipping thread-budget test"
        end
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
