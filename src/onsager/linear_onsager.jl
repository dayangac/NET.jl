"""
Module for linear Onsager relations.
"""
module LinearOnsager

using Plots

# Manual logging implementation
function log_info(message::String)
    println("[INFO] $message")
end

# Manual matrix-vector multiplication
function matvec_mult(A::Array{Float64, 2}, x::Array{Float64, 1})::Array{Float64, 1}
    m, n = size(A)
    result = zeros(Float64, m)
    for i in 1:m
        for j in 1:n
            result[i] += A[i, j] * x[j]
        end
    end
    return result
end

"""
    struct OnsagerMatrix

Represents a phenomenological coefficient matrix for linear Onsager relations.
- L: 2D array for a uniform matrix or 3D array for a spatially varying matrix.
"""
mutable struct OnsagerMatrix
    L::Union{Array{Float64, 2}, Array{Float64, 3}}
end

"""
    compute_fluxes(L::OnsagerMatrix, F::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) based on forces (F) using the phenomenological coefficient matrix (L).
- L: An `OnsagerMatrix` type representing the phenomenological coefficients.
- F: 2D array representing forces at each spatial point.
Returns: A 2D array of computed fluxes.
"""
function compute_fluxes(L::OnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions(L, F)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)
    for k in 1:num_points
        if ndims(L.L) == 2
            J[:, k] = matvec_mult(L.L, F[:, k])
        else
            J[:, k] = matvec_mult(L.L[:, :, k], F[:, k])
        end
    end
    log_info("Flux computation complete for $num_points points.")
    return J
end



"""
    validate_dimensions(L::OnsagerMatrix, F::Array{Float64, 2})

Ensures that the input dimensions of the matrix L and the force array F match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions(L::OnsagerMatrix, F::Array{Float64, 2})
    num_vars, num_points = size(F)
    if ndims(L.L) == 2
        size(L.L, 1) == num_vars || throw(DimensionMismatch("Number of variables in L and F must match."))
    elseif ndims(L.L) == 3
        (size(L.L, 1) == num_vars && size(L.L, 3) == num_points) ||
            throw(DimensionMismatch("L must match the number of variables and points in F."))
    else
        throw(ArgumentError("L must be a 2D or 3D array."))
    end
    return nothing
end

"""
    visualize_fluxes(J::Array{Float64, 2}, title_name::String)

Creates a heatmap visualization of the fluxes for easy interpretation.
- J: 2D array of fluxes at each spatial point.
"""
function visualize_fluxes(J::Array{Float64, 2}, title_name::String)
    heatmap(J, xlabel="Grid Points", ylabel="Flux Variables",
            title=title_name, color=:viridis, colorbar=true)
end

export OnsagerMatrix, compute_fluxes, validate_dimensions, visualize_fluxes

end # module LinearOnsager