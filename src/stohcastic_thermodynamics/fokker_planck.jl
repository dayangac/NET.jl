"""
Module for Fokker-Planck equation solutions.
"""
module FokkerPlanck

"""
    solve_fokker_planck(x_range, t_range, D, drift)

Solves the Fokker-Planck equation for a given drift field and diffusion coefficient.

# Arguments
- x_range::Vector{Float64}: Spatial range for the solution.
- t_range::Vector{Float64}: Time range for the solution.
- D::Float64: Diffusion coefficient.
- drift::Function: Drift function, `drift(x)`.

# Returns
- solution::Matrix{Float64}: Probability distribution over space and time.
"""
function solve_fokker_planck(x_range, t_range, D, drift)
    # Calculate step sizes manually
    dx = length(x_range) > 1 ? x_range[2] - x_range[1] : 1.0
    dt = length(t_range) > 1 ? t_range[2] - t_range[1] : 1.0


    P = zeros(length(x_range), length(t_range))
    P[:, 1] .= exp.(-x_range.^2 / 2) ./ sqrt(2π)  


    for t in 1:(length(t_range) - 1)
        for i in 2:(length(x_range) - 1)

            drift_term = -drift(x_range[i]) * (P[i+1, t] - P[i-1, t]) / (2 * dx)
            diff_term = D * (P[i+1, t] - 2 * P[i, t] + P[i-1, t]) / dx^2
            P[i, t+1] = P[i, t] + dt * (drift_term + diff_term)
        end
    end

    return P
end

export solve_fokker_planck

end # module FokkerPlanck
