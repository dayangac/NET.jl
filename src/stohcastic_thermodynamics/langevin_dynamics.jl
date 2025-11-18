"""
Module for Langevin dynamics simulations.
"""
module LangevinDynamics

# Physical constants - k_B is defined in the main NET module
# Access k_B from parent module
k_B = 1.380649e-23  # Boltzmann constant (J/K)

# Manual random number generator (Linear Congruential Generator)
const RNG_A = 1664525
const RNG_C = 1013904223
const RNG_M = 2^32
rng_state = UInt64(12345)  # Seed

function randn()::Float64
    global rng_state
    # Box-Muller transform for normal distribution
    u1 = rand()
    u2 = rand()
    z = sqrt(-2.0 * log(u1 + 1e-10)) * cos(2.0 * π * u2)
    return z
end

function rand()::Float64
    global rng_state
    rng_state = (RNG_A * rng_state + RNG_C) % RNG_M
    return Float64(rng_state) / RNG_M
end

function randn(n::Int)::Array{Float64, 1}
    result = zeros(Float64, n)
    for i in 1:n
        result[i] = randn()
    end
    return result
end

"""
    langevin_dynamics(num_steps, dt, gamma, T, mass, potential, grad_potential; dims=1)

Simulates Langevin dynamics for a particle in a potential field with thermal noise, damping, and deterministic forces.

# Arguments
- `num_steps::Int`: Number of simulation steps.
- `dt::Float64`: Time step for the simulation (seconds).
- `gamma::Float64`: Damping coefficient (kg/s).
- `T::Float64`: Temperature (Kelvin).
- `mass::Float64`: Mass of the particle (kg).
- `potential::Function`: Function defining the potential energy field.
- `grad_potential::Function`: Function defining the gradient (force) of the potential.
- `dims::Int`: Dimensionality of the simulation (default: 1).

# Returns
- `time::Vector{Float64}`: Time vector.
- `positions::Matrix{Float64}`: Position of the particle at each time step for each dimension.
- `velocities::Matrix{Float64}`: Velocity of the particle at each time step for each dimension.
"""
function langevin_dynamics(num_steps::Int, dt::Float64, gamma::Float64, T::Float64, mass::Float64,
                           potential::Function, grad_potential::Function; dims::Int = 1)

    time = collect(0:dt:(num_steps - 1) * dt)
    positions = zeros(Float64, num_steps, dims)
    velocities = zeros(Float64, num_steps, dims)

    sigma = sqrt(2 * k_B * T * gamma / mass)
    friction_factor = exp(-gamma * dt / mass)

    positions[1, :] .= randn(dims)
    velocities[1, :] .= sigma * randn(dims)

    for i in 2:num_steps
        x = positions[i - 1, :]
        v = velocities[i - 1, :]

        random_noise = sigma * sqrt(1 - friction_factor^2) * randn(dims)
        v_half = v .* friction_factor .+ 0.5 * dt .* (-grad_potential(x) ./ mass) .+ random_noise
        x_new = x .+ dt .* v_half
        v_new = v_half .+ 0.5 * dt .* (-grad_potential(x_new) ./ mass)

        positions[i, :] .= x_new
        velocities[i, :] .= v_new
    end

    return time, positions, velocities
end

export langevin_dynamics

end # module LangevinDynamics
