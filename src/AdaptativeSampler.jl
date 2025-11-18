"""
    adaptive_sampling(f, interval, budget; initial_points=[], tolerance=1e-3)

Performs adaptive sampling of a function `f(x)` over the specified `interval` (a tuple `(a, b)`).
The function uses a budget constraint (maximum number of samples) and refines sampling adaptively
based on the function's behavior. Returns a vector of sampled points and their corresponding values.

# Arguments
- `f`: The function to sample.
- `interval`: A tuple `(a, b)` specifying the sampling interval.
- `budget`: The maximum number of samples allowed.

# Keyword Arguments
- `initial_points`: A vector of initial x-values to start sampling (default: `[]`).
- `tolerance`: The tolerance for detecting significant changes in `f(x)` (default: `1e-3`).

# Returns
A tuple `(x_samples, y_samples)` where:
- `x_samples`: A vector of sampled x-values.
- `y_samples`: A vector of corresponding f(x)-values.
"""
function adaptive_sampling(f, interval::Tuple{Number, Number}, budget::Int; initial_points=[], tolerance=1e-3)
    a, b = interval

    function compute_midpoint(x1, x2)
        return (x1 + x2) / 2
    end
    # Start with initial points or the interval endpoints
    if isempty(initial_points)
        x_samples = [a, compute_midpoint(a, b), b]
    else
        x_samples = initial_points[a .<=  initial_points .<= b]
        push!(x_samples, a)
        push!(x_samples, b) 
        sort!(x_samples)
        unique!(x_samples) 
    end
    if length(x_samples) < 2
        error("At least two points are required for adaptive sampling.")
    end
    if length(x_samples) == 2
        # If only two points, sample the midpoint
        mid_x = compute_midpoint(x_samples[1], x_samples[2])
        insert(x_samples, 1, mid_x)
    end
    y_samples = [f(x) for x in x_samples]
    remaining_budget = budget - length(x_samples)

    while remaining_budget > 0
            errors = map(1:length(x_samples)-2) do idx
                y1 = y_samples[idx]
                y2 = y_samples[idx+2]
                y_mid = y_samples[idx+1]
                y_estimate = y1 + (y2 - y1) / (x_samples[idx+2] - x_samples[idx]) * (x_samples[idx+1] - x_samples[idx]) 
                return abs(y_mid - y_estimate)
            end
            idx_max_error = argmax(errors)
            if errors[idx_max_error] > tolerance
                if x_samples[idx_max_error+1] - x_samples[idx_max_error] > x_samples[idx_max_error+2] - x_samples[idx_max_error+1]
                    new_x = compute_midpoint(x_samples[idx_max_error], x_samples[idx_max_error+1])
                    insert_position = idx_max_error+1
                else
                    new_x = compute_midpoint(x_samples[idx_max_error+1], x_samples[idx_max_error+2])
                    insert_position = idx_max_error+2
                end
                insert!(x_samples, insert_position, new_x)
                insert!(y_samples, insert_position, f(new_x))
            end
        remaining_budget -= 1
    end

    return x_samples, y_samples
end
