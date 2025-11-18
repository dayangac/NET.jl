"""
Manual ODE solver implementation to replace DifferentialEquations dependency.
Implements Runge-Kutta 4th order method for solving ordinary differential equations.
"""
module ODESolver

"""
    rk4_solve(f, u0, tspan, dt; reltol=1e-8, abstol=1e-8)

Solves an ODE using 4th order Runge-Kutta method.

# Arguments
- `f`: Function f(u, p, t) that returns du/dt
- `u0`: Initial condition
- `tspan`: Tuple (t_start, t_end)
- `dt`: Time step size
- `reltol`: Relative tolerance (optional)
- `abstol`: Absolute tolerance (optional)

# Returns
- `t`: Time points
- `u`: Solution values at each time point
"""
function rk4_solve(f, u0, tspan, dt; reltol=1e-8, abstol=1e-8)
    t_start, t_end = tspan
    t = collect(t_start:dt:t_end)
    if t[end] != t_end
        push!(t, t_end)
    end
    
    n = length(t)
    u = zeros(typeof(u0), n)
    u[1] = u0
    
    for i in 1:(n-1)
        dt_actual = t[i+1] - t[i]
        k1 = f(u[i], nothing, t[i])
        k2 = f(u[i] + dt_actual * k1 / 2, nothing, t[i] + dt_actual / 2)
        k3 = f(u[i] + dt_actual * k2 / 2, nothing, t[i] + dt_actual / 2)
        k4 = f(u[i] + dt_actual * k3, nothing, t[i] + dt_actual)
        
        u[i+1] = u[i] + (dt_actual / 6) * (k1 + 2*k2 + 2*k3 + k4)
    end
    
    return (t, u)
end

"""
    ODESolution

Simple struct to hold ODE solution results, compatible with DifferentialEquations interface.
"""
struct ODESolution
    t::Vector{Float64}
    u::Vector{Float64}
    
    function ODESolution(t, u)
        new(t, u)
    end
end

"""
    solve_ode(f, u0, tspan; dt=0.01, reltol=1e-8, abstol=1e-8)

High-level interface for solving ODEs, similar to DifferentialEquations.solve.

# Arguments
- `f`: Function f(u, p, t) that returns du/dt
- `u0`: Initial condition
- `tspan`: Tuple (t_start, t_end)
- `dt`: Time step size (default: 0.01)
- `reltol`: Relative tolerance (optional)
- `abstol`: Absolute tolerance (optional)

# Returns
- `ODESolution`: Solution object with t and u fields
"""
function solve_ode(f, u0, tspan; dt=0.01, reltol=1e-8, abstol=1e-8)
    t, u = rk4_solve(f, u0, tspan, dt; reltol=reltol, abstol=abstol)
    return ODESolution(t, u)
end

export rk4_solve, solve_ode, ODESolution

end # module ODESolver

