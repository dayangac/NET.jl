# NET.jl

[![Build Status](https://github.com/dayangac/NET.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dayangac/NET.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Docs

A Julia library for Non-Equilibrium Thermodynamics (NET) analysis and simulation.

## Table of Contents

- [Introduction](#introduction)
- [Onsager Relations](#onsager-relations)
  - [Linear Onsager Relations](#linear-onsager-relations)
  - [Non-Linear Onsager Relations](#non-linear-onsager-relations)
  - [Stochastic Onsager Relations](#stochastic-onsager-relations)
  - [Time-Dependent Onsager Relations](#time-dependent-onsager-relations)
- [Entropy Production](#entropy-production)
  - [Local Entropy Production](#local-entropy-production)
  - [Time-Dependent Entropy Production](#time-dependent-entropy-production)
  - [Entropy Production Rate](#entropy-production-rate)
  - [Entropy Optimization](#entropy-optimization)
- [Stochastic Thermodynamics](#stochastic-thermodynamics)
  - [Stochastic Driven Systems](#stochastic-driven-systems)
  - [Fluctuation Theorems](#fluctuation-theorems)
  - [Langevin Dynamics](#langevin-dynamics)
  - [Fokker-Planck Equation](#fokker-planck-equation)
  - [Path Integrals](#path-integrals)
  - [Non-Equilibrium Potential](#non-equilibrium-potential)
  - [Noise-Induced Bifurcations](#noise-induced-bifurcations)
- [Geometrothermodynamics](#geometrothermodynamics)
  - [Contact Geometry](#contact-geometry)
  - [Legendre Transformations](#legendre-transformations)

---

## Introduction

NET.jl is a comprehensive Julia library designed for working with Non-Equilibrium Thermodynamics principles. The library provides tools for visualization, simulation, and computation of non-equilibrium thermodynamic systems. It encompasses three primary domains:

1. **Onsager Relations** - Describing flux-force relationships in thermodynamic systems
2. **Entropy Production** - Quantifying energy dissipation and irreversibility
3. **Stochastic Thermodynamics** - Analyzing systems influenced by random fluctuations

---

## Onsager Relations

Onsager relations express the fundamental relationship between thermodynamic fluxes and forces in non-equilibrium systems. Forces such as temperature gradients and concentration differences drive fluxes, which determine system behaviors including heat transfer and particle diffusion.

The library provides four implementations: linear, non-linear, time-dependent, and stochastic Onsager relations.

### Linear Onsager Relations

Linear Onsager relations describe systems where fluxes are linearly proportional to thermodynamic forces. The implementation supports both 2D and 3D Onsager matrices.

#### Structure Definition

```julia
mutable struct OnsagerMatrix
    L::Union{Array{Float64, 2}, Array{Float64, 3}}
end
```

The dimensionality determines the system type:
- **2D matrices**: Simple linear Onsager relationships
- **3D matrices**: Spatial Onsager relations accounting for local thermal properties

#### Flux Computation

The fundamental relation is given by:

$$\mathbf{J} = \mathbf{L} \cdot \mathbf{F}$$

where $\mathbf{J}$ is the flux vector, $\mathbf{L}$ is the Onsager matrix, and $\mathbf{F}$ is the force vector.

```julia
function compute_fluxes(L::OnsagerMatrix, F::Array{Float64, 2})
    validate_dimensions(L, F)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)
    
    for k in 1:num_points
        if ndims(L.L) == 2
            J[:, k] = L.L * F[:, k]
        else
            J[:, k] = L.L[:, :, k] * F[:, k]
        end
    end
    
    log_info("Flux computation complete for $num_points points.")
    return J
end
```

For 3D matrices, each spatial point possesses individual thermal properties represented as a sub-matrix $\mathbf{L}_{:,:,k}$ at position $k$.

### Non-Linear Onsager Relations

Non-linear Onsager relations arise in systems far from equilibrium where higher-order coupling terms become significant. These systems require spatial resolution of thermal properties.

#### Structure Definition

```julia
mutable struct NonLinearOnsagerMatrix
    L::Array{Float64, 3}
end

mutable struct NonLinearFluxSystem
    matrix::NonLinearOnsagerMatrix
    linear_forces::Array{Float64, 2}
    nonlinear_forces::Vector{Array{Float64, 2}}
end
```

#### Total Flux Computation

The total flux includes both linear and non-linear contributions:

$$\mathbf{J}_{\text{total}} = \mathbf{L} \cdot \mathbf{F}_{\text{linear}} + \sum_{i} \mathbf{L} \cdot \mathbf{F}_{\text{nonlinear},i}$$

```julia
function compute_total_fluxes(system::NonLinearFluxSystem)
    L = system.matrix.L
    F = system.linear_forces
    F_nl_list = system.nonlinear_forces

    validate_dimensions_multinonlinear(system)

    num_vars, num_points = size(F)
    J_linear = zeros(num_vars, num_points)
    J_total = zeros(num_vars, num_points)

    # Compute linear contribution
    for k in 1:num_points
        J_linear[:, k] = L[:, :, k] * F[:, k] 
    end
    J_total .= J_linear 

    # Add non-linear contributions
    for F_nl in F_nl_list  
        for k in 1:num_points
            J_total[:, k] .+= L[:, :, k] * F_nl[:, k]
        end
    end

    return J_total
end
```

### Stochastic Onsager Relations

Stochastic Onsager relations extend spatial relations to include explicit noise terms, accounting for thermal fluctuations and random perturbations in the system.

#### Structure Definition

```julia
mutable struct StochasticOnsagerMatrix
    L::Array{Float64, 3}
end
```

#### Stochastic Flux Computation

The stochastic flux incorporates both deterministic and random components:

$$\mathbf{J} = \mathbf{L} \cdot \mathbf{F} + \boldsymbol{\xi}$$

where $\boldsymbol{\xi}$ represents the noise term.

```julia
function compute_stochastic_fluxes(L::StochasticOnsagerMatrix, 
                                   F::Array{Float64, 2}, 
                                   noise::Array{Float64, 2})
    validate_dimensions_stochastic(L, F, noise)
    num_vars, num_points = size(F)
    J = zeros(num_vars, num_points)

    for k in 1:num_points
        J[:, k] = L.L[:, :, k] * F[:, k] + noise[:, k]
    end

    log_info("Stochastic flux computation complete for $num_points points with external noise.")
    return J
end
```

### Time-Dependent Onsager Relations

Time-dependent Onsager relations describe systems where transport coefficients vary with time rather than space, applicable to systems with temporal gradients.

#### Temporal Flux Computation

```julia
function compute_time_dependent_fluxes(L::TimeDependentOnsagerMatrix, 
                                       F::Array{Float64, 2})
    validate_dimensions_time_dependent(L, F)
    num_vars, num_times = size(F)
    J = zeros(num_vars, num_times)

    for t in 1:num_times
        J[:, t] = L.L[:, :, t] * F[:, t]
    end

    log_info("Time-dependent flux computation complete for $num_times time points.")
    return J
end
```

**Note**: The 3D matrix structure allows 2D matrices to be multiplied with 3D matrices when the first two dimensions match. Each depth slice is treated as an independent 2D matrix operation.

---

## Entropy Production

Entropy production quantifies the rate at which a system generates entropy, representing irreversible energy dissipation. It is driven by thermodynamic forces and provides insight into the degree of irreversibility in non-equilibrium processes.

### Local Entropy Production

Local entropy production is computed as the inner product of flux and force vectors:

$$\sigma = \mathbf{J} \cdot \mathbf{X}$$

where $\sigma$ is the local entropy production rate, $\mathbf{J}$ is the flux vector, and $\mathbf{X}$ is the thermodynamic force vector.

```julia
function compute_local_entropy_production(J::Array{Float64, 2}, X::Array{Float64, 2})
    validate_dimensions_entropy(J, X)
    entropy_production = sum(J .* X, dims=1)  
    log_info("Local entropy production calculation complete for $(size(J, 2)) grid points.")
    return vec(entropy_production)  
end
```

#### Dimension Validation

```julia
function validate_dimensions_entropy(J::Array{Float64, 2}, X::Array{Float64, 2})
    if size(J) != size(X)
        throw(DimensionMismatch("J and X must have the same dimensions."))
    end
end
```

#### Visualization

```julia
function visualize_local_entropy_production_heatmap(σ::Array{Float64, 1}, 
                                                     grid_points::Vector, 
                                                     title_name::String)
    σ_matrix = reshape(σ, 1, length(σ))
    
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
```

### Time-Dependent Entropy Production

Time-dependent entropy production analyzes how entropy generation evolves dynamically. The system is modeled using ordinary differential equations.

The rate of entropy production is given by:

$$\frac{d\sigma}{dt} = \mathbf{J}(t) \cdot \mathbf{X}(t)$$

```julia
function analyze_time_dependent_entropy_production(J_func::Function, 
                                                    X_func::Function, 
                                                    t_span::Tuple{Float64, Float64})
    function dσ_dt(p, t)
        J = p * J_func(t)
        X = X_func(t)
        return sum(J .* X)
    end

    σ_0 = 0.0  
    prob = ODEProblem(dσ_dt, σ_0, t_span)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    return sol
end
```

The function accepts time-dependent flux and force functions, along with the time span for integration. The parameter `p` allows for scaling coefficients to be applied to the flux.

### Entropy Production Rate

The entropy production rate quantifies the change in entropy production over a finite time interval:

$$\frac{\Delta \sigma}{\Delta t} = \frac{\sigma(t_2) - \sigma(t_1)}{\Delta t}$$

```julia
function compute_entropy_rate(J1::Array{Float64, 2}, X1::Array{Float64, 2},
                              J2::Array{Float64, 2}, X2::Array{Float64, 2}, 
                              Δt::Float64)
    validate_dimensions_entropy(J1, X1)
    validate_dimensions_entropy(J2, X2)

    entropy_initial = sum(J1 .* X1, dims=1)
    entropy_final = sum(J2 .* X2, dims=1)

    Δentropy = sum(entropy_final) - sum(entropy_initial)
    entropy_rate = Δentropy / Δt

    log_info("Entropy rate calculation complete: rate = $(entropy_rate) per unit time.")
    return entropy_rate
end
```

### Entropy Optimization

Entropy optimization employs gradient descent to minimize system entropy, useful for finding minimum entropy production states.

#### Gradient Computation

The gradients with respect to flux and force are:

$$\nabla_J \sigma = \mathbf{X}, \quad \nabla_X \sigma = \mathbf{J}$$

```julia
function entropy_gradient(J::Array{Float64, 2}, X::Array{Float64, 2})
    ∇J = X  
    ∇X = J  
    return (∇J, ∇X)
end

function gradient_descent(J::Array{Float64, 2}, X::Array{Float64, 2}, 
                         learning_rate::Float64, max_iterations::Int)
    for i in 1:max_iterations
        ∇J, ∇X = entropy_gradient(J, X)
        J -= learning_rate * ∇J
        X -= learning_rate * ∇X
        
        @show compute_entropy(J, X)
    end
    return (J, X)
end
```

---

## Stochastic Thermodynamics

Stochastic thermodynamics extends classical thermodynamics to systems influenced by random fluctuations. It provides a framework for analyzing energy, entropy, and noise in small-scale systems governed by stochastic dynamics.

### Stochastic Driven Systems

This module simulates stochastic environments with physically realistic conditions, including:

- **Diffusion process**: Random motion due to thermal fluctuations
- **Drift**: Systematic motion caused by dominant forces
- **Stochastic growth**: Random growth and decay in particle properties
- **Non-linear noise dynamics**: Higher-order interactions affecting system behavior

Two numerical methods are implemented for solving stochastic differential equations (SDEs):

#### Euler-Maruyama Method

The Euler-Maruyama scheme is suitable for systems with weak noise:

$$X_{t+\Delta t} = X_t + f(X_t)\Delta t + g(X_t)\sqrt{\Delta t}\,\xi_t$$

where:
- $f(X_t)$ is the drift function
- $g(X_t)$ is the noise amplitude
- $\xi_t \sim \mathcal{N}(0,1)$ is Gaussian white noise (Wiener process)
- $\Delta t$ is the time step

```julia
if method == "Euler-Maruyama"
    trajectories[t, r] = trajectories[t - 1, r] +
                         drift(trajectories[t - 1, r]) * dt +
                         noise_amplitude(trajectories[t - 1, r]) * sqrt(dt) * ξ
```

This method provides first-order strong convergence and is computationally efficient for small noise amplitudes.

#### Milstein Method

The Milstein scheme includes second-order noise corrections for improved accuracy in non-linear systems:

$$X_{t+\Delta t} = X_t + f(X_t)\Delta t + g(X_t)\sqrt{\Delta t}\,\xi_t + \frac{1}{2}g(X_t)g'(X_t)\left(\xi_t^2 - 1\right)\Delta t$$

where $g'(X_t)$ is the derivative of the noise amplitude function.

```julia
elseif method == "Milstein"
    g = noise_amplitude(trajectories[t - 1, r])
    g_prime = (noise_amplitude(trajectories[t - 1, r] + 1e-6) - g) / 1e-6 
    trajectories[t, r] = trajectories[t - 1, r] +
                         drift(trajectories[t - 1, r]) * dt +
                         g * sqrt(dt) * ξ +
                         0.5 * g * g_prime * (ξ^2 - 1) * dt
```

The Milstein method captures non-linear noise interactions and provides higher-order accuracy, derived using Itô calculus.

### Fluctuation Theorems

Fluctuation theorems describe statistical properties of thermodynamic fluctuations in non-equilibrium systems. The library implements two fundamental theorems:

#### Jarzynski Equality

The Jarzynski equality relates non-equilibrium work to equilibrium free energy:

$$\left\langle e^{-\beta W} \right\rangle = e^{-\beta \Delta F}$$

where:
- $\beta = 1/(k_B T)$ is the inverse temperature
- $W$ is the work performed on the system
- $\Delta F$ is the equilibrium free energy difference
- $\langle \cdot \rangle$ denotes ensemble average

```julia
function jarzynski_equality(work_values::Vector{Float64}, T::Float64)
    beta = 1 / (k_B * T)

    # Log-sum-exp trick for numerical stability
    max_work = maximum(-beta * work_values)
    avg_exp_work = exp(-max_work) * mean(exp.(-beta * work_values .+ max_work))
    
    if avg_exp_work <= 0 || isnan(avg_exp_work)
        println("Warning: avg_exp_work is non-positive or NaN; adjusting to avoid -Inf result.")
        avg_exp_work = 1e-15
    end
    
    ΔF_estimated = -log(avg_exp_work) / beta
    return ΔF_estimated
end
```

The implementation uses the log-sum-exp trick (RealSoftMax) to maintain numerical stability when computing exponential averages.

#### Crooks Fluctuation Theorem

Crooks' theorem relates forward and reverse work distributions:

$$\frac{P_F(W)}{P_R(-W)} = e^{\beta(W - \Delta F)}$$

where:
- $P_F(W)$ is the probability of observing work $W$ in the forward process
- $P_R(-W)$ is the probability of observing work $-W$ in the reverse process
- $\Delta F$ is the free energy difference

The theorem enables estimation of equilibrium free energy from non-equilibrium measurements at the crossing point where $P_F(W) = P_R(-W)$.

```julia
function crooks_theorem(work_forward::Vector{Float64}, 
                       work_reverse::Vector{Float64}, T::Float64)
    beta = 1 / (k_B * T)  

    log_ratios = Float64[]  
  
    for (w_f, w_r) in zip(work_forward, work_reverse)
        log_numerator = beta * w_f
        log_denominator = beta * w_r

        if isfinite(log_numerator) && isfinite(log_denominator) && log_denominator > -700
            push!(log_ratios, log_numerator - log_denominator)
        end
    end

    if isempty(log_ratios)
        println("Warning: No valid log-ratios found; returning NaN.")
        return NaN
    end

    mean_log_ratio = mean(log_ratios)
    mean_ratio = exp(mean_log_ratio)
    return mean_ratio
end
```

The function returns the probability ratio of forward to reverse work, quantifying the relative likelihood of each process.

### Langevin Dynamics

Langevin dynamics models particle motion in a thermal bath, combining deterministic forces, friction, and stochastic thermal noise. The Langevin equation is:

$$m\frac{dv}{dt} = -\nabla U(x) - \gamma v + \sqrt{2\gamma k_B T}\,\xi(t)$$

where:
- $m$ is the particle mass
- $v$ is velocity
- $U(x)$ is the potential energy
- $\gamma$ is the friction coefficient
- $k_B$ is Boltzmann's constant
- $T$ is temperature
- $\xi(t)$ is white noise with $\langle\xi(t)\xi(t')\rangle = \delta(t-t')$

The thermal noise satisfies the fluctuation-dissipation theorem, ensuring that energy dissipated through friction is balanced by thermal fluctuations, maintaining equilibrium at temperature $T$.

#### BAOAB Integrator

The BAOAB integrator is a splitting scheme that preserves equilibrium properties while accurately propagating Langevin dynamics. The method decomposes each timestep into five operations:

1. **B**: Half velocity update using forces
2. **A**: Position update using velocity
3. **O**: Stochastic momentum update (Ornstein-Uhlenbeck process)
4. **A**: Second position update
5. **B**: Second half velocity update

```julia
function langevin_dynamics(num_steps::Int, dt::Float64, gamma::Float64, T::Float64, 
                          mass::Float64, potential::Function, grad_potential::Function; 
                          dims::Int = 1)

    time = collect(0:dt:(num_steps - 1) * dt)
    positions = zeros(Float64, num_steps, dims)
    velocities = zeros(Float64, num_steps, dims)

    sigma = sqrt(2 * k_B * T * gamma / mass)
    friction_factor = exp(-gamma * dt / mass)

    # Initialize
    positions[1, :] .= randn(dims)
    velocities[1, :] .= sigma * randn(dims)

    for i in 2:num_steps
        x = positions[i - 1, :]
        v = velocities[i - 1, :]

        # BAOAB scheme
        random_noise = sigma * sqrt(1 - friction_factor^2) * randn(dims)
        v_half = v .* friction_factor .+ 0.5 * dt .* (-grad_potential(x) ./ mass) .+ random_noise
        x_new = x .+ dt .* v_half
        v_new = v_half .+ 0.5 * dt .* (-grad_potential(x_new) ./ mass)

        positions[i, :] .= x_new
        velocities[i, :] .= v_new
    end

    return time, positions, velocities
end
```

The BAOAB method ensures accurate sampling of the Boltzmann distribution while maintaining stability and geometric properties of the dynamics.

### Fokker-Planck Equation

The Fokker-Planck equation describes the evolution of the probability density function (PDF) for stochastic systems:

$$\frac{\partial P(x,t)}{\partial t} = -\frac{\partial}{\partial x}\left[f(x)P(x,t)\right] + D\frac{\partial^2 P(x,t)}{\partial x^2}$$

where:
- $P(x,t)$ is the probability density at position $x$ and time $t$
- $f(x)$ is the drift function
- $D$ is the diffusion coefficient

The first term represents deterministic drift, while the second term accounts for diffusive spreading.

```julia
function solve_fokker_planck(x_range, t_range, D, drift)
    dx = step(x_range)
    dt = step(t_range)

    P = zeros(length(x_range), length(t_range))
    P[:, 1] .= exp.(-x_range.^2 / 2) ./ sqrt(2π)  # Initial Gaussian distribution

    for t in 1:(length(t_range) - 1)
        for i in 2:(length(x_range) - 1)
            # Drift term: -∂[f(x)P]/∂x
            drift_term = -drift(x_range[i]) * (P[i+1, t] - P[i-1, t]) / (2 * dx)
            
            # Diffusion term: D∂²P/∂x²
            diff_term = D * (P[i+1, t] - 2 * P[i, t] + P[i-1, t]) / dx^2
            
            P[i, t+1] = P[i, t] + dt * (drift_term + diff_term)
        end
    end

    return P
end
```

The implementation uses finite difference methods with forward Euler time integration.

### Path Integrals

Path integrals provide a framework for analyzing the most probable trajectories in stochastic systems. The action functional quantifies path likelihood:

$$S[x(t)] = \int_0^T \left[\frac{1}{2D}\left(\frac{dx}{dt} - f(x)\right)^2\right] dt$$

where:
- $S[x(t)]$ is the action for path $x(t)$
- $f(x)$ is the drift
- $D$ is the diffusion coefficient
- $T$ is the time duration

Paths minimizing the action correspond to the most probable trajectories. The path probability is proportional to $e^{-S[x(t)]}$.

```julia
function calculate_action(path, drift, diffusion)
    action = 0.0
    for i in 1:(length(path) - 1)
        dx = path[i+1] - path[i]
        drift_term = drift(path[i]) * dx / diffusion
        action += (dx^2 / (2 * diffusion)) + drift_term
    end
    return action
end
```

The function discretizes the integral, computing the action for a given trajectory through the system's phase space.

### Non-Equilibrium Potential

The non-equilibrium potential provides a landscape representation of system dynamics, revealing stable states, basins of attraction, and transition pathways:

$$\Phi(x) = -\int^x \frac{f(x')}{D} dx'$$

where:
- $\Phi(x)$ is the non-equilibrium potential
- $f(x')$ is the drift at position $x'$
- $D$ is the diffusion coefficient

Minima of $\Phi(x)$ correspond to stable states, while maxima represent unstable transition states.

```julia
function calculate_non_equilibrium_potential(grid, drift, diffusion)
    potential = zeros(length(grid))
    for i in 2:length(grid)
        dx = grid[i] - grid[i-1]
        potential[i] = potential[i-1] + drift(grid[i]) * dx / diffusion
    end
    return potential
end
```

### Noise-Induced Bifurcations

Bifurcations represent qualitative transitions in system dynamics. Noise-induced bifurcations occur when stochastic fluctuations alter stability properties or create new attractors.

The deterministic system has fixed points where $f(x) = 0$. When noise is introduced, the effective dynamics change. Critical transitions occur when:

$$|f(x)| \lesssim \sigma_{\text{noise}}$$

where the deterministic drift becomes comparable to noise intensity.

Additional stability analysis considers the effective potential in the presence of noise:

$$\Phi_{\text{eff}}(x) = -\int \frac{f(x)}{D} dx - \frac{1}{2}\log g(x)$$

where $g(x)$ is the noise amplitude function (if state-dependent).

```julia
function detect_bifurcation(drift, noise_level, x_range)
    bifurcation_points = []
    for x in x_range
        if abs(drift(x)) < noise_level
            push!(bifurcation_points, x)
        end
    end
    return bifurcation_points
end
```

The function identifies regions where noise can induce qualitative changes in system behavior, marking potential bifurcation points.

---

## Geometrothermodynamics

Geometrothermodynamics applies differential geometry to thermodynamic systems, providing a geometric framework for understanding phase transitions and equilibrium properties.

### Contact Geometry

Contact geometry extends symplectic geometry to describe non-equilibrium thermodynamic phase spaces. A contact manifold is characterized by a contact form that encodes relationships between thermodynamic potentials and variables.

#### Contact Form

The contact form $\Theta$ on a thermodynamic phase space is defined as:

$$\Theta = dE - \sum_{\alpha} I^{\alpha} dX^{\alpha}$$

where:
- $E$ is the thermodynamic potential (e.g., internal energy)
- $X^{\alpha}$ are extensive variables (e.g., volume, entropy)
- $I^{\alpha}$ are intensive variables (e.g., pressure, temperature)

```julia
function contact_form(E::Function, X::Vector{Symbol}, I::Vector{Symbol})
    n = length(X)
    dE = Symbol("∂$E/∂$(join(X, ','))") 
    dX = ["d$(X[i])" for i in 1:n]
    Θ = string(dE, " - ", join([string(I[i], " * ", dX[i]) for i in 1:n], " - "))
    return Θ
end
```

#### Non-Degeneracy Condition

A contact manifold must satisfy the non-degeneracy condition:

$$\Theta \wedge (d\Theta)^n \neq 0$$

where $n$ is the dimensionality of the system. This ensures the manifold has proper geometric structure.

```julia
function non_degeneracy_condition(contact_form::String, n::Int)
    return "$(contact_form) ∧ (d$(contact_form))^$n ≠ 0"
end
```

#### Thermodynamic Metric

The thermodynamic metric defines the geometric structure of the manifold:

$$g = g_{\alpha\beta} dx^{\alpha} dx^{\beta}$$

```julia
function thermodynamic_metric(g_ab::Matrix, variables::Vector{Symbol})
    dx = ["d$var" for var in variables]
    metric = sum([g_ab[i, j] * "($dx[$i])($dx[$j])" 
                  for i in 1:length(dx), j in 1:length(dx)])
    return metric
end
```

#### Entropy-Based Metric

For systems described by entropy $S(U, X^{\alpha})$, the metric incorporates entropy derivatives:

$$g = -\frac{\partial^2 S}{\partial U^2} dU^2 + \sum_{\alpha,\beta} \frac{\partial^2 S}{\partial X^{\alpha} \partial X^{\beta}} dX^{\alpha} dX^{\beta}$$

```julia
function entropy_based_metric(S::Function, U::Symbol, X::Vector{Symbol})
    d2S_dU2 = Symbol("∂²$S/∂$U²")
    d2S_dXiXj = [Symbol("∂²$S/∂$(X[i])∂$(X[j])") 
                 for i in 1:length(X), j in 1:length(X)]
    
    metric_U = "-($d2S_dU2 * d$U^2)"
    metric_X = join([string(d2S_dXiXj[i], " d", X[i], " d", X[j]) 
                     for i in 1:length(X)], " + ")
    
    return metric_U * " + " * metric_X
end
```

This metric captures the curvature of the thermodynamic phase space and relates geometric properties to physical observables.

### Legendre Transformations

Legendre transformations are mathematical operations that convert between different thermodynamic potentials by changing the set of independent variables. These transformations preserve physical content while providing different perspectives on the system.

For a potential $E(X^{\alpha})$, the Legendre transform with respect to variable $X^{k}$ yields:

$$\tilde{E}(I^{1}, \ldots, I^{k}, X^{k+1}, \ldots, X^{n}) = E - I^{k}X^{k}$$

where $I^{k} = \partial E / \partial X^{k}$ is the conjugate intensive variable.

Common thermodynamic potentials related by Legendre transforms include:
- Internal energy $U(S,V,N)$
- Helmholtz free energy $F(T,V,N) = U - TS$
- Enthalpy $H(S,P,N) = U + PV$
- Gibbs free energy $G(T,P,N) = U - TS + PV$

These transformations are fundamental to thermodynamic formalism and enable analysis from different experimental conditions (e.g., constant temperature vs. constant entropy).

---

## Installation

```julia
using Pkg
Pkg.add("NET")
```

## Usage

```julia
using NET

# Example: Linear Onsager relation
L = OnsagerMatrix([1.0 0.5; 0.5 1.0])
F = [1.0 2.0; 0.5 1.5]
J = compute_fluxes(L, F)

# Example: Langevin dynamics
time, pos, vel = langevin_dynamics(
    1000,                              # number of steps
    0.01,                              # time step
    1.0,                               # friction coefficient
    300.0,                             # temperature (K)
    1.0,                               # mass
    x -> 0.5 * x^2,                   # harmonic potential
    x -> x;                           # gradient
    dims=1
)
```

## Contributing

Contributions are welcome. Please submit pull requests or open issues for bugs and feature requests.

## License

This project is licensed under the MIT License.

## References

For theoretical background on non-equilibrium thermodynamics, stochastic processes, and geometrothermodynamics, refer to standard texts in statistical mechanics and differential geometry.

---

**Author**: Ed
**Institution**: Hisar School IdeaLab FabLab, Istanbul, Turkey
**Date**: November 2024

## License

MIT License - see [LICENSE](LICENSE) file for details.
