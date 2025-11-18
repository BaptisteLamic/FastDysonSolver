using CpuId
function select_BLAS()
    # Check if the CPU supports AVX2
    if cpuvendor() == :Intel
        println("Intel CPU detected")
        @eval using MKL
    elseif  Sys.isapple()
        println("Apple CPU detected")
        @eval using AppleAccelerate
    end
end