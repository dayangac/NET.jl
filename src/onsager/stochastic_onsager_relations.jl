"""
Module for stochastic Onsager relations.
"""
module StochasticOnsagerRelations

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
    struct StochasticOnsagerMatrix

Represents a phenomenological matrix for stochastic Onsager relations.
- L: 3D array where L[:,:,k] is the matrix of coefficients at grid point or time step k.
"""
mutable struct StochasticOnsagerMatrix
    L::Array{Float64, 3}
end

"""
    compute_stochastic_fluxes(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2}) -> Array{Float64, 2}

Computes the fluxes (J) with stochastic contributions at each point based on the forces (F) using the phenomenological matrix (L).
- L: A `StochasticOnsagerMatrix` type representing the stochastic phenomenological coefficients.
- F: 2D array representing forces at each point.
- noise: 2D array representing noise contributions at each point.
Returns: A 2D array of computed fluxes with stochastic contributions.
"""
function compute_stochastic_fluxes(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})
    validate_dimensions_stochastic(L, F, noise)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)

    for k in 1:num_points
        J[:, k] = matvec_mult(L.L[:, :, k], F[:, k]) + noise[:, k]
    end

    log_info("Stochastic flux computation complete for $num_points points with external noise.")
    return J
end

"""
    validate_dimensions_stochastic(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})

Ensures that the input dimensions of the matrix L, the force array F, and the noise array match.
Throws an error if the dimensions are not compatible.
"""
function validate_dimensions_stochastic(L::StochasticOnsagerMatrix, F::Array{Float64, 2}, noise::Array{Float64, 2})
    num_vars, num_points = size(F)
    if size(L.L, 1) != num_vars || size(L.L, 3) != num_points
        throw(DimensionMismatch("L must match the number of variables and points in F."))
    end
    if size(noise) != (num_vars, num_points)
        throw(DimensionMismatch("Noise dimensions must match the dimensions of F (variables × points)."))
    end
end

"""
    visualize_stochastic_fluxes(J::Array{Float64, 2}, points::Vector, title_name::String)

Creates a heatmap to visualize stochastic fluxes over points (e.g., space or time).
- J: 2D array of fluxes at each point.
- points: Vector representing the x-axis (e.g., spatial or time points).
"""
function visualize_stochastic_fluxes(J::Array{Float64, 2}, points::Vector, title_name::String)
    p = heatmap(
        points,
        1:size(J, 1),
        J,
        xlabel="Grid Points",
        ylabel="Flux Variables",
        title=title_name,
        color=:inferno,
        colorbar=true
    )
    display(p)
end

export StochasticOnsagerMatrix, compute_stochastic_fluxes, 
       validate_dimensions_stochastic, visualize_stochastic_fluxes

end # module StochasticOnsagerRelations
