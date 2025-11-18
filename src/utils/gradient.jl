"""
Manual gradient computation to replace ForwardDiff dependency.
Implements finite difference method for computing gradients.
"""
module Gradient

"""
    gradient(f, x; h=1e-6)

Computes the gradient of a function f at point x using finite differences.

# Arguments
- `f`: Function that takes a vector and returns a scalar
- `x`: Point at which to compute the gradient
- `h`: Step size for finite differences (default: 1e-6)

# Returns
- `grad`: Gradient vector
"""
function gradient(f, x; h=1e-6)
    n = length(x)
    grad = zeros(Float64, n)
    fx = f(x)
    
    for i in 1:n
        x_plus = copy(x)
        x_plus[i] += h
        fx_plus = f(x_plus)
        grad[i] = (fx_plus - fx) / h
    end
    
    return grad
end

export gradient

end # module Gradient

