# This interpolates on intervall [-1,1]×[-1,1]
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

# This interpolates on intervall [-1,1]×[-1,1]×[-1,1]
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

