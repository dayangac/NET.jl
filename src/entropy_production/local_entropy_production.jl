"""
Module for local entropy production calculations.
"""
module LocalEntropyProduction

using Plots

# Manual logging implementation
function log_info(message::String)
    println("[INFO] $message")
end

# Manual linear algebra functions
function sum_dims(A::Array{Float64, 2}, dim::Int)
    """Manual implementation of sum(A, dims=dim)"""
    if dim == 1
        # Sum along rows (first dimension), result is a row vector
        m, n = size(A)
        result = zeros(Float64, 1, n)
        for j in 1:n
            for i in 1:m
                result[1, j] += A[i, j]
            end
        end
        return result
    elseif dim == 2
        # Sum along columns (second dimension), result is a column vector
        m, n = size(A)
        result = zeros(Float64, m, 1)
        for i in 1:m
            for j in 1:n
                result[i, 1] += A[i, j]
            end
        end
        return result
    else
        throw(ArgumentError("dim must be 1 or 2"))
    end
end

function vec_manual(A::Array{Float64, 2})
    """Manual implementation of vec() - converts matrix to vector column-wise"""
    m, n = size(A)
    result = zeros(Float64, m * n)
    idx = 1
    for j in 1:n
        for i in 1:m
            result[idx] = A[i, j]
            idx += 1
        end
    end
    return result
end

function reshape_manual(A::Array{Float64, 1}, dims::Tuple{Int, Int})
    """Manual implementation of reshape() - reshapes vector to matrix"""
    m, n = dims
    if length(A) != m * n
        throw(DimensionMismatch("Cannot reshape array of length $(length(A)) to size ($m, $n)"))
    end
    result = zeros(Float64, m, n)
    idx = 1
    for j in 1:n
        for i in 1:m
            result[i, j] = A[idx]
            idx += 1
        end
    end
    return result
end

function elementwise_multiply_sum(A::Array{Float64, 2}, B::Array{Float64, 2})
    """Manual implementation of sum(A .* B)"""
    m, n = size(A)
    if size(B) != (m, n)
        throw(DimensionMismatch("A and B must have the same size"))
    end
    result = 0.0
    for i in 1:m
        for j in 1:n
            result += A[i, j] * B[i, j]
        end
    end
    return result
end

"""
    compute_local_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2}) -> Array{Float64, 1}

Calculates the local entropy production at each grid point based on the fluxes (J) and corresponding thermodynamic forces (X).
- J: 2D array of fluxes at each point, where each column represents a flux vector at a grid point.
- X: 2D array of thermodynamic forces at each point, where each column represents a force vector at a grid point.
Returns: A 1D array of entropy production values for each grid point.
"""
function compute_local_entropy_prozduction(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)
    # Manual element-wise multiplication and sum along rows
    m, n = size(J)
    entropy_production_2d = zeros(Float64, 1, n)
    for j in 1:n
        for i in 1:m
            entropy_production_2d[1, j] += J[i, j] * X[i, j]
        end
    end
    log_info("Local entropy production calculation complete for $(size(J, 2)) grid points.")
    return vec_manual(entropy_production_2d)  
end

"""
    validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})

Ensures that the input dimensions of the fluxes (J) and forces (X) match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    if size(J) != size(X)
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
    return nothing
end

"""
    visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, grid_points::Vector, title_name::String)

Creates a heatmap to visualize local entropy production across grid points.
- σ: 1D array of entropy production values at each grid point.
- grid_points: Vector representing the grid points for the x-axis.
"""
function visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, grid_points::Vector, title_name::String)

    σ_matrix = reshape_manual(σ, (1, length(σ)))
    
    p = heatmap(
        grid_points,
        1:1,  
        σ_matrix,
        xlabel="Grid Points",
        ylabel="Entropy Production",
        title=title_name,
        color=:plasma,
        colorbar=true
    )
    display(p) 
end

export compute_local_entropy_prozduction, validate_dimensions_entropy,
       visualize_local_entropy_production_heatmap

end # module LocalEntropyProduction
