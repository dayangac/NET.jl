"""
Module for entropy production rate calculations.
"""
module EntropyProductionRate

# Manual logging implementation
function log_info(message::String)
    println("[INFO] $message")
end

# Manual sum function for 2D arrays along dimension
function sum_dims(A::Array{Float64, 2}, dim::Int)
    """Manual implementation of sum(A, dims=dim)"""
    if dim == 1
        m, n = size(A)
        result = zeros(Float64, 1, n)
        for j in 1:n
            for i in 1:m
                result[1, j] += A[i, j]
            end
        end
        return result
    elseif dim == 2
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

# Manual sum function for arrays (handles both 1D and 2D)
function sum_manual(A::Array{Float64})
    """Manual implementation of sum() - works for both 1D and 2D arrays"""
    result = 0.0
    if ndims(A) == 1
        for x in A
            result += x
        end
    elseif ndims(A) == 2
        m, n = size(A)
        for i in 1:m
            for j in 1:n
                result += A[i, j]
            end
        end
    else
        throw(ArgumentError("sum_manual only supports 1D and 2D arrays"))
    end
    return result
end

"""
    compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                         J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64) -> Float64

Calculates the rate of irreversible entropy production by comparing entropy production at two different states over time Δt.
- J1, X1: 2D arrays of fluxes and thermodynamic forces at the initial state.
- J2, X2: 2D arrays of fluxes and thermodynamic forces at the final state.
- Δt: Time interval between the two states.
Returns: The rate of entropy production in the system per unit time.
"""
function compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                              J2::Array{Float64, 2}, X2::Array{Float64, 2}, Δt::Float64)
    validate_dimensions_entropy(J1, X1)
    validate_dimensions_entropy(J2, X2)

    # Calculate local entropy production at initial and final states
    m1, n1 = size(J1)
    entropy_initial_2d = zeros(Float64, 1, n1)
    for j in 1:n1
        for i in 1:m1
            entropy_initial_2d[1, j] += J1[i, j] * X1[i, j]
        end
    end
    
    m2, n2 = size(J2)
    entropy_final_2d = zeros(Float64, 1, n2)
    for j in 1:n2
        for i in 1:m2
            entropy_final_2d[1, j] += J2[i, j] * X2[i, j]
        end
    end

    # Calculate the change in entropy production and normalize by time
    Δentropy = sum_manual(entropy_final_2d) - sum_manual(entropy_initial_2d)
    entropy_rate = Δentropy / Δt

    log_info("Entropy rate calculation complete: rate = $(entropy_rate) per unit time.")
    return entropy_rate
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

export compute_entropy_rate
# Note: validate_dimensions_entropy is not exported to avoid ambiguity with LocalEntropyProduction

end # module EntropyProductionRate
