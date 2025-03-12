"""
    gen_cont_data_3D(n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0, σ3::Float64=5.0; pos_only::Bool=true, mode::String="abs", min_val::Float64=1.0e-10, is_periodic::Bool=false, mean_zero::Bool=false)

Generate a 3D array of continuous data with optional periodicity, positivity, and Gaussian filtering.

# Arguments
- `n_elem::Int=100`: Number of elements along each dimension of the output array.
- `σ1::Float64=5.0`: Standard deviation for Gaussian filter in the first dimension.
- `σ2::Float64=5.0`: Standard deviation for Gaussian filter in the second dimension.
- `σ3::Float64=5.0`: Standard deviation for Gaussian filter in the third dimension.
- `pos_only::Bool=true`: If `true`, ensures the output array has only positive values.
- `mode::String="abs"`: Method to ensure positivity:
  - `"abs"`: Take the absolute value of the random data before filtering.
  - `"exp"`: Apply exponential function to the filtered data.
- `min_val::Float64=1.0e-10`: Minimum value added to ensure positivity.
- `is_periodic::Bool=false`: If `true`, generates periodic data using circular padding.
- `mean_zero::Bool=false`: If `true`, adjusts the data to have zero mean (only if `pos_only` is `false`).

# Returns
- A 3D array of size `(n_elem, n_elem, n_elem)` with the generated continuous data, scaled so that its mean is 1.

# Notes
- The output is scaled such that its mean is 1, making it suitable for further computations.
- If `is_periodic` is `false`, the function generates a larger array (3 times the size in the first two dimensions) and extracts the central `(n_elem, n_elem, n_elem)` portion to avoid boundary effects.
- Requires the `Statistics` package for mean calculations and assumes `gauß_filter` is defined elsewhere when `is_periodic=true`.
- The function uses `imfilter` for Gaussian filtering in the non-periodic case, which is a Julia-native function.
"""

function gen_cont_data_3D( n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0, σ3::Float64=5.0,; pos_only::Bool=true, mode::String="abs", min_val::Float64=1.0e-10, is_periodic::Bool=false, mean_zero::Bool=false, use_scipy::Bool=true)
    # Initialize array with random entries
    if is_periodic
        A = randn(n_elem, n_elem, n_elem)
        if pos_only
            if mode == "abs"
                A = abs.(A)

                A = gauß_filter(A, (σ1, σ2, σ3), "circular", use_scipy)
            elseif mode == "exp"
                A = gauß_filter(A, (σ1, σ2, σ3), "circular", use_scipy)
                A = exp.(A)
            else
                A = gauß_filter(A, (σ1, σ2, σ3), "circular", use_scipy)
            end
            A .+=  -minimum(A)+ min_val
        else
            A = gauß_filter(A, (σ1, σ2, σ3), "circular", use_scipy)
            if mean_zero
                A .-= Statistics.mean(A)
            end
        end

    
    else # The nonperiodic case
        A = randn(3*n_elem, 3*n_elem, n_elem)    
        if pos_only
            if mode == "abs"
                A = abs.(A)
                A = imfilter(A, Kernel.gaussian((σ1, σ2, σ3)), "replicate")
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            elseif mode == "exp"
                A = imfilter(A, Kernel.gaussian((σ1, σ2, σ3)), "replicate")
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem,n_elem+1:2*n_elem]
                A = exp.(A)
            else
                A = imfilter(A, Kernel.gaussian((σ1, σ2, σ3)), "replicate")
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            end
            A .+= - minimum(A)+ min_val
        else # In case negative values are allowed
            A = imfilter(A, Kernel.gaussian((σ1, σ2, σ3)), "replicate")
            A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            if mean_zero
                A .-= Statistics.mean(A)
            end
        end
    end

    #Scale it to useable value
    mean_A = Statistics.mean(A)
    return (1 / mean_A) * A
end