"""
    interpolate_array_2D(arr::Array{Float64, 2}, positive::Bool=false)

Create an interpolation function for a 2D array `arr` over the interval [-1,1]×[-1,1] or [0,1]×[0,1].

# Arguments
- `arr::Array{Float64, 2}`: The input 2D array to interpolate. Must be square (n × n).
- `positive::Bool=false`: If `true`, interpolate over [0,1]×[0,1]; otherwise, over [-1,1]×[-1,1].

# Returns
- A function that takes a 2D point `x` (a vector with two elements) and returns the interpolated value from `arr`.

# Notes
- The interpolation is linear and uses `Interpolations.interpolate` with `Gridded(Linear())`.
- The input array must be square; otherwise, an assertion error is raised.
- When `positive=true`, the interpolation maps [0,1]×[0,1] to the array indices.
- When `positive=false`, it maps [-1,1]×[-1,1] to the array indices.
"""
function interpolate_array_2D(arr::Array{Float64, 2}, positive::Bool=false)
    # Ensure the input array is n x n
    @assert size(arr, 1) == size(arr, 2) "The input array must be square (n x n)."

    # Define the range of the original array indices
    n = size(arr, 1)
    xs = 1:n
    ys = 1:n
    
    # Create an interpolation object
    itp = Interpolations.interpolate((xs, ys), arr, Gridded(Linear()))
    
    # Define a function to map the interval [-1, 1] to the array index range [1, n]
    if positive # If true interpolates on intervall [0,1]×[0,1]
        return x -> itp(1 + x[1] * (n - 1), 1 + x[2] * (n - 1))
    else
        return x -> itp(1 + (0.5*x[1]+ 0.5) * (n - 1), 1 + (0.5*x[2]+0.5) * (n - 1))
    end
end
"""
    interpolate_array_3D(arr::Array{Float64, 3}, positive::Bool=false)

Create an interpolation function for a 3D array `arr` over the interval [-1,1]×[-1,1]×[-1,1] or [0,1]×[0,1]×[0,1].

# Arguments
- `arr::Array{Float64, 3}`: The input 3D array to interpolate. Must be cubic (n × n × n).
- `positive::Bool=false`: If `true`, interpolate over [0,1]×[0,1]×[0,1]; otherwise, over [-1,1]×[-1,1]×[-1,1].

# Returns
- A function that takes a 3D point `x` (a vector with three elements) and returns the interpolated value from `arr`.

# Notes
- The interpolation is linear and uses `Interpolations.interpolate` with `Gridded(Linear())`.
- The input array must be cubic; otherwise, an assertion error is raised.
- When `positive=true`, the interpolation maps [0,1]×[0,1]×[0,1] to the array indices.
- When `positive=false`, it maps [-1,1]×[-1,1]×[-1,1] to the array indices.
"""
function interpolate_array_3D(arr::Array{Float64, 3}, positive::Bool=false)
    # Ensure the input array is n x n x n
    @assert size(arr, 1) == size(arr, 2) == size(arr, 3) "The input array must be cubic (n x n x n)."
    # Define the range of the original array indices
    n = size(arr, 1)
    xs = 1:n
    ys = 1:n
    zs = 1:n    # Create an interpolation object
    itp = Interpolations.interpolate((xs, ys, zs), arr, Gridded(Linear()))      
    # Define a function to map the interval [-1, 1] to the array index range [1, n]
    if positive # If true interpolates on intervall [0,1]×[0,1]×[0,1]
        return x -> itp(1 + x[1] * (n - 1), 1 + x[2] * (n - 1), 1 + x[3] * (n - 1))
    else
        return x -> itp(1 + (0.5*x[1]+ 0.5) * (n - 1), 1 + (0.5*x[2]+0.5) * (n - 1), 1 + (0.5*x[3]+0.5) * (n - 1))
    end
end

