export gauß_filter
    # This function uses white noise and filters it.
include("Domain2D.jl")
export gen_cont_data_2D, gen_discrete_data_2D, combine_cont_dicrete_2D

    # Scipy's filter seems to be better. Images.jl has this annoying gridlike anomalies:
"""
    gauß_filter(A, σ_g, mode::String="circular", scipy::Bool=false)

Apply a Gaussian filter to the input array `A` with standard deviations `σ_g`.

# Arguments
- `A`: The input array to be filtered.
- `σ_g`: A tuple or array specifying the standard deviations for the Gaussian filter along each axis.
- `mode::String="circular"`: The mode for handling boundaries. Options are:
  - `"circular"`: Wrap-around using circular padding.
  - `"replicate"`: Replicate the edge values.
- `scipy::Bool=false`: If `true`, use SciPy's `gaussian_filter` for filtering; otherwise, use Julia's `imfilter`.

# Returns
- The filtered array.

# Notes
- If `scipy=true`, the function requires the SciPy library to be installed and accessible via Julia's PyCall or similar interfaces.
- The `mode` parameter is translated to SciPy's mode when `scipy=true`:
  - `"circular"` corresponds to `"wrap"`.
  - `"replicate"` corresponds to `"nearest"`.
- When `scipy=false`, the function uses Julia's `imfilter` with the specified `mode`.
"""
function gauß_filter(A,σ_g, mode::String="circular", scipy::Bool=false)
    if scipy == true
        if mode == "circular"
            scipy_mode = "wrap"
        elseif mode =="replicate"
            scipy_mode = "nearest"
        end
        return SciPy.ndimage.gaussian_filter(A, σ_g, mode =scipy_mode)
    else
        return imfilter(A, Kernel.gaussian(σ_g), mode)
    end
end


    


include("Boundary2D.jl") # Not really necessary
export gen_cont_data_1D, gen_single_points_1D
include("Domain3D.jl")
export gen_cont_data_3D
include("Interpolate.jl")
export interpolate_array_2D, interpolate_array_3D