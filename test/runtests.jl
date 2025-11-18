using DrWatson, Test
@quickactivate "FastDysonSolver"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
includet("./../src/Junction.jl")
println("Starting tests")
ti = time()

@testset "FastDysonSolver tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
