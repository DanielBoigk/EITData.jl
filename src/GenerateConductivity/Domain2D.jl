
"""
    gen_cont_data_2D(n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0; pos_only::Bool=true, mode::String="abs", min_val::Float64=1.0e-10, is_periodic::Bool=false, mean_zero::Bool=false, use_scipy::Bool=true)

Generate a 2D array of continuous data with optional periodicity, positivity, and Gaussian filtering.

# Arguments
- `n_elem::Int=100`: Number of elements along each dimension of the output array.
- `σ1::Float64=5.0`: Standard deviation for Gaussian filter in the first dimension.
- `σ2::Float64=5.0`: Standard deviation for Gaussian filter in the second dimension.
- `pos_only::Bool=true`: If true, ensures the output array has only positive values.
- `mode::String="abs"`: Method to ensure positivity:
  - `"abs"`: Take the absolute value of the random data before filtering.
  - `"exp"`: Apply exponential function to the filtered data.
- `min_val::Float64=1.0e-10`: Minimum value added to ensure positivity.
- `is_periodic::Bool=false`: If true, generates periodic data using circular padding.
- `mean_zero::Bool=false`: If true, adjusts the data to have zero mean (only if `pos_only` is false).
- `use_scipy::Bool=true`: If true, uses SciPy for Gaussian filtering via the `gauß_filter` function.

# Returns
- A 2D array of size `(n_elem, n_elem)` with the generated continuous data.

# Notes
- The output is scaled such that its mean (or mean of absolute values if `pos_only` is false) is approximately 1.
- If `is_periodic` is false, the function generates a larger array (3 times the size) and extracts the central `(n_elem, n_elem)` portion to avoid boundary effects.
- Requires the `Statistics` package for mean calculations and assumes `gauß_filter` is defined elsewhere.
"""

function gen_cont_data_2D( n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0,; pos_only::Bool=true, mode::String="abs", min_val::Float64=1.0e-10, is_periodic::Bool=false, mean_zero::Bool=false, use_scipy::Bool=true)
    # Initialize array with random entries
    if is_periodic
        A = randn(n_elem, n_elem)
        if pos_only
            if mode == "abs"
                A = abs.(A)
                A = gauß_filter(A, (σ1, σ2), "circular", use_scipy)
            elseif mode == "exp"
                A = gauß_filter(A, (σ1, σ2), "circular", use_scipy)
                A = exp.(A)
            else
              
                A = gauß_filter(A, (σ1, σ2), "circular", use_scipy)
            end
            A .+=  -minimum(A)+ min_val
        else
            A = gauß_filter(A, (σ1, σ2), "circular", use_scipy)
            if mean_zero
                A .-= Statistics.mean(A)
            end
        end

    
    else # The nonperiodic case
        A = randn(3*n_elem, 3*n_elem)    
        if pos_only
            if mode == "abs"
                A = abs.(A)
                #A = imfilter(A, Kernel.gaussian((σ1, σ2)), "replicate")
                A = gauß_filter(A, (σ1, σ2), "replicate", use_scipy)
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            elseif mode == "exp"
                #A = imfilter(A, Kernel.gaussian((σ1, σ2)), "replicate")
                A = gauß_filter(A, (σ1, σ2), "replicate", use_scipy)
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem]
                A = exp.(A)
            else
                #A = imfilter(A, Kernel.gaussian((σ1, σ2)), "replicate")
                A = gauß_filter(A, (σ1, σ2), "replicate", use_scipy)
                A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            end
            A .+= - minimum(A)+ min_val
        else # In case negative values are allowed
            #A = imfilter(A, Kernel.gaussian((σ1, σ2)), "replicate")
            A = gauß_filter(A, (σ1, σ2), "replicate", use_scipy)
            A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem]
            if mean_zero
                A .-= Statistics.mean(A)
            end
        end
    end
    if pos_only
        #Scale it to useable value
        mean_A = Statistics.mean(A)
        
    else
        mean_A = Statistics.mean(abs.(A))
    end
    return (1 / mean_A) * A
end



"""
    gen_discrete_data_2D(n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0; pos_only::Bool=true, mode::String="abs", threshold::Float64=0.0, is_periodic::Bool=false, set_one_zero::Bool=false, use_scipy::Bool=true)

Generate a 2D array with discrete values based on a threshold applied to filtered random data.

# Arguments
- `n_elem::Int=100`: Number of elements along each dimension of the output array.
- `σ1::Float64=5.0`: Standard deviation for Gaussian filter in the first dimension.
- `σ2::Float64=5.0`: Standard deviation for Gaussian filter in the second dimension.
- `pos_only::Bool=true`: If true, ensures the discrete values are positive.
- `mode::String="abs"`: Not used in this function (included for compatibility with `gen_cont_data_2D`).
- `threshold::Float64=0.0`: Threshold value to determine discrete levels.
- `is_periodic::Bool=false`: If true, generates periodic data using circular padding.
- `set_one_zero::Bool=false`: If true, forces one of the discrete values to be zero.
- `use_scipy::Bool=true`: If true, uses SciPy for Gaussian filtering via the `gauß_filter` function.

# Returns
- A 2D array of size `(n_elem, n_elem)` with discrete values (either `a` or `b`).

# Notes
- The discrete values `a` and `b` are randomly generated, with `a` set to 0 if `set_one_zero` is true.
- If `pos_only` is true, both `a` and `b` are non-negative.
- If `is_periodic` is false, the function generates a larger array (3 times the size) and extracts the central `(n_elem, n_elem)` portion to avoid boundary effects.
- Assumes `gauß_filter` is defined elsewhere.
"""
function gen_discrete_data_2D(n_elem::Int=100, σ1::Float64=5.0, σ2::Float64=5.0,; pos_only::Bool=true, mode::String="abs", threshold::Float64=0.0, is_periodic::Bool=false, set_one_zero::Bool=false, use_scipy::Bool=true)
if set_one_zero
    a = 0
    b = randn()
else
    a = randn()
    b = randn()
end
if pos_only
    a = abs(a)
    b = abs(b)
end    
if set_one_zero
    a = 0
end
# Initialize array with random entries
if is_periodic
    A = randn(n_elem, n_elem)

    A = gauß_filter(A, (σ1, σ2), "circular", use_scipy)
else
    A = randn(3*n_elem, 3*n_elem)

    A = gauß_filter(A, (σ1, σ2), "replicate", use_scipy)
    A = A[n_elem+1:2*n_elem,n_elem+1:2*n_elem]
end
for i in eachindex(A)
    A[i] = A[i] > threshold ? a : b
end

return A
end


"""
    combine_cont_dicrete_2D(D::Array{Float64,2}, C1::Array{Float64,2}, C2::Array{Float64,2}, is_zero::Bool=false, ensure_positive::Bool=false)

Combine discrete and continuous data arrays based on the values in the discrete array.

# Arguments
- `D::Array{Float64,2}`: The discrete data array.
- `C1::Array{Float64,2}`: The first continuous data array to add.
- `C2::Array{Float64,2}`: The second continuous data array to add.
- `is_zero::Bool=false`: If true, uses zero as the threshold for combining.
- `ensure_positive::Bool=false`: Not used in this function (included for potential compatibility).

# Returns
- A 2D array resulting from adding `C1` or `C2` to `D` based on the values in `D`.

# Notes
- If `is_zero` is true, adds `C1` where `D` is zero and `C2` elsewhere.
- If `is_zero` is false, adds `C1` where `D` equals its minimum value and `C2` elsewhere.
- The function assumes that `D`, `C1`, and `C2` have the same dimensions.
"""

function combine_cont_dicrete_2D(D::Array{Float64,2},C1::Array{Float64,2},C2::Array{Float64,2}, is_zero::Bool = false, ensure_positive = false)
    A = copy(D)
    values = unique(A)
    values = sort(values)
    if is_zero
        for i in eachindex(A)
            A[i] = A[i] == 0.0 ? A[i]+C1[i] : A[i]+C2[i]
        end
    else
        for i in eachindex(A)
            A[i] = A[i] == values[1] ? A[i]+C1[i] : A[i]+C2[i]
        end
    end
    return A
end