function select_BLAS()
    # Optional BLAS backends: used when installed in the active environment,
    # skipped otherwise. MKL has no Apple-Silicon build and AppleAccelerate is
    # macOS-only, so neither is a hard dependency of the project.
    if Sys.isapple() && Base.identify_package("AppleAccelerate") !== nothing
        println("Apple CPU detected, using AppleAccelerate")
        @eval using AppleAccelerate
    elseif Base.identify_package("MKL") !== nothing
        println("MKL detected in the environment, using MKL")
        @eval using MKL
    else
        println("Using default OpenBLAS")
    end
end
