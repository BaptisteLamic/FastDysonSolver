using BenchmarkTools

using LinearAlgebra
using HssMatrices
using NonEquilibriumGreenFunction
using StaticArrays

using CairoMakie
using Colors
using ColorSchemes
using LaTeXStrings

using SpecialFunctions
using Symbolics
using SymbolicNumericIntegration
using Interpolations
using DSP

using QuadGK

#First we define a structure and its contructors to keep together all the simulation parameters
struct JunctionParameters
    δt::Float64 #timestep
    T::Float64  #simulation length
    Γl::Float64 #left tunneling rate
    Γr::Float64 #right tunneling rate
    β::Float64  #inverse temperature
    Δ::Float64
    ϕl # left phase
    ϕr # right phase
    τ::Float64 # Dot apodization time
    ε_τ::Float64 # Apodization tolerance
    n::Int64 # Apodization order
end

function similar(p::JunctionParameters; kwargs...)
    # Get all field names and current values
    fields = fieldnames(typeof(p))
    values = [getfield(p, f) for f in fields]
    # Override with any keyword arguments provided
    for (k, v) in kwargs
        idx = findfirst(isequal(k), fields)
        if idx !== nothing
            values[idx] = v
        else
            error("Unknown field name: $k")
        end
    end
    return JunctionParameters(; zip(fields, values)...)
end
function JunctionParameters(; δt, T, Γl=1, Γr=1, β=1000, Δ=1.0, ϕl = x->0, ϕr = x->0, τ = 40T, ε_τ = 0.9, n = 0) # by default there is no apodization
    return JunctionParameters(δt, T, Γl, Γr, β, Δ, ϕl, ϕr, τ, ε_τ, n)
end

#The leafsize should be adapated to the problem at hand. 
#If it is too large, or too small, the solver is slowdown. 
default_compression() = HssCompression( atol = 1E-8, rtol = 1E-8,leafsize = 64, kest = 20)

axis(p::JunctionParameters) = 0:p.δt:p.T
N(p::JunctionParameters) = length(axis(p))
function N_from_observation(observation)
    N(observation[:param])
end
σ0() = @SMatrix [1.0 0.0; 0.0 1.0]
σx() = @SMatrix [0 1.0; 1.0 0.0]
σz() = @SMatrix [1.0 0.0; 0.0 -1.0]

function get_apodization_function(p::JunctionParameters)
    # Returns the apodization function for the given parameters.
    # The order n is used to define the sharpness of the cutoff, and the time constant is set by ε_τ and τ.
    if p.n == 0
        return (t, tp) -> 1.0
    end
    return (t, tp) -> exp(-((abs(t-tp)/_evaluate_apodization_time_constante(p.ε_τ,p.τ,p.n))^(2p.n)))
end
function _evaluate_apodization_time_constante(ε_τ,τ,n)
    return τ * log(1/ε_τ)^(-1/2n)
end

function compute_retarded_lead_green_function(p::JunctionParameters; cpr)
    g_R_lead_delta = discretize_dirac(axis(p), t -> -1im* pi * σ0(), compression=cpr)
    g_R_lead_continuous = discretize_retardedkernel(axis(p),
        (t, tp) ->  pi*(p.Δ * besselj0(p.Δ*(t-tp)) * σx()+ 1im * p.Δ * besselj1(p.Δ*(t-tp)) * σ0()),
        compression=cpr, stationary=true)
    g_lead = g_R_lead_delta + g_R_lead_continuous
    return g_lead
end

 
function compute_GR(p::JunctionParameters;cpr)
    g_R_lead = compute_retarded_lead_green_function(p::JunctionParameters; cpr=cpr)
    coupling_left = discretize_dirac(axis(p), t ->  sqrt(p.Γl/pi) * exp(1im * σz() * p.ϕl(t) / 2) * σz(), compression=cpr)
    coupling_right = discretize_dirac(axis(p), t -> sqrt(p.Γr/pi) * exp(1im* σz() * p.ϕr(t) / 2) * σz(), compression=cpr)
    Σ_R_left = coupling_left' * g_R_lead * coupling_left
    Σ_R_right = coupling_right' * g_R_lead * coupling_right
    Σ_R = Σ_R_left + Σ_R_right
    apodization_function = get_apodization_function(p)
    g = discretize_retardedkernel(axis(p), (t, tp) -> ComplexF64(-1im)*exp(-1im * σz() * 0 * (t-tp)) * apodization_function(t, tp), compression=cpr, stationary=true)
    G_R = solve_dyson(g, g * Σ_R)
    return (; g_R_lead,  g, G_R, Σ_R_left, Σ_R_right, Σ_R, coupling_left, coupling_right)
end

function simulate_junction(p::JunctionParameters; cpr = default_compression() )
    results_GR = compute_GR(p,cpr = cpr)
    g_R_lead = results_GR[:g_R_lead]
    coupling_left = results_GR[:coupling_left]
    coupling_right = results_GR[:coupling_right]
    G_R = results_GR[:G_R]
    ρ = discretize_acausalkernel(axis(p), (t, tp) -> thermal_kernel(t - tp, p.β) * σ0() .|> ComplexF64,
    stationary=true, compression=cpr)
    g_lead_kinetic = g_R_lead *  ρ  - ρ * g_R_lead'

    Σ_K_left =  compress!(coupling_left' * g_lead_kinetic * coupling_left)
    Σ_K_right = compress!(coupling_right' * g_lead_kinetic * coupling_right)
    Σ_K = Σ_K_left + Σ_K_right
    G_K = compress!(G_R * Σ_K * G_R') #+ dot_correction
    return (;results_GR..., G_K, Σ_K_left, Σ_K_right)
end

function simulate_thermal_equilibrium_junction(p::JunctionParameters; cpr = default_compression() )
    results_GR = compute_GR(p,cpr = cpr)
    g_R_lead = results_GR[:g_R_lead]
    coupling_left = results_GR[:coupling_left]
    coupling_right = results_GR[:coupling_right]
    G_R = results_GR[:G_R]
    ρ = discretize_acausalkernel(axis(p), (t, tp) -> thermal_kernel(t - tp, p.β) * σ0() .|> ComplexF64,
    stationary=true, compression=cpr)
    g_lead_kinetic = g_R_lead *  ρ  - ρ * g_R_lead'
    Σ_K_left =  compress!(coupling_left' * g_lead_kinetic * coupling_left)
    Σ_K_right = compress!(coupling_right' * g_lead_kinetic * coupling_right)
    G_K = G_R *  ρ  - ρ * G_R'
    return (;results_GR..., G_K, Σ_K_left, Σ_K_right)
end



function compute_average_current(results)
    #First we build the expression
    function f(G_R,G_K,Σ_R,Σ_K)
        _Σ_R_time_G_R = Threads.@spawn Σ_R*G_R
        _G_R_time_Σ_R = Threads.@spawn G_R*Σ_R
        _G_R_time_Σ_K = Threads.@spawn G_R*Σ_K  
        _G_K_time_adjoint_Σ_R = Threads.@spawn G_K*adjoint(Σ_R)
        _Σ_K_time_adjoint_G_R = Threads.@spawn Σ_K*adjoint(G_R)
        _Σ_R_time_G_K = Threads.@spawn Σ_R*G_K

        Σ_R_time_G_R = fetch(_Σ_R_time_G_R)
        G_R_time_Σ_R = fetch(_G_R_time_Σ_R)
        G_R_time_Σ_K = fetch(_G_R_time_Σ_K)
        G_K_time_adjoint_Σ_R = fetch(_G_K_time_adjoint_Σ_R)
        Σ_K_time_adjoint_G_R = fetch(_Σ_K_time_adjoint_G_R)
        Σ_R_time_G_K = fetch(_Σ_R_time_G_K)
        
        return -((1//2)*((adjoint(Σ_R_time_G_R)) - Σ_R_time_G_R) + G_R_time_Σ_K + G_K_time_adjoint_Σ_R + -(Σ_K_time_adjoint_G_R + Σ_R_time_G_K) + (3//2)*(G_R_time_Σ_R  - adjoint(G_R_time_Σ_R )))
    end
    I_avr_op = f(results[:G_R], results[:G_K], results[:Σ_R_left], results[:Σ_K_left])
    #We have to take the trace on the Keldysh space by hand. 
    dg = diag(matrix(I_avr_op))
    (dg[1:2:end] .- dg[2:2:end])
end

function getAverage_JosepshonFrequency(p::JunctionParameters)
    #Warning, assume constant phase direction
    return abs((p.ϕl(p.T) -  p.ϕr(p.T))/(2π*p.T)) 
end

function get_ωj_filter(p)
    n_TJ = 2/getAverage_JosepshonFrequency(p)
    n_steps = round(Int, n_TJ / p.δt)
    kernel = [ 1 / n_steps for i in 1:n_steps]
   return function(x)
        data = DSP.conv(x, kernel)
        return data[1:end-length(kernel)+1]
   end
end

function extract_dc(I,p)
    ωj_filter = get_ωj_filter(p)
    filtred = ωj_filter(I)
    N = length(filtred)
    return mean(filtred[div(N,2):end])
end
 
function single_run(p,cpr = default_compression())
    result = simulate_junction(p,cpr = cpr)
    return result
end

function single_run_current(p,cpr = default_compression())
    result = simulate_junction(p,cpr = cpr)
    current = compute_average_current(result)
    return (;p, I = real.(current))
end

function single_run_thermal_equilibrium_current(p,phi,cpr; Tmax = nothing)
    T = 3 * p.β
    if Tmax isa Number
        T = min(T,Tmax)
        if T == Tmax
            @warn "T/p.β small"
        end
    end
    _p = similar(p, ϕl=t -> -phi/2, ϕr=t -> phi/2)
    result = simulate_thermal_equilibrium_junction(_p,cpr = cpr)
    current = compute_average_current(result)
    Ith = real.(current)
    return (;p = _p, I = Ith)
end


function measure_single_run(p,cpr)
    return @benchmark single_run($p,$cpr)
end
