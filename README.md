# NET.jl

[![Build Status](https://github.com/dayangac/NET.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dayangac/NET.jl/actions/workflows/CI.yml?query=branch%3Amain)

**NET.jl** - Non-Equilibrium Thermodynamics Library for Julia

A comprehensive Julia library for working with Non-Equilibrium Thermodynamics (NET) principles, including Onsager Relations, Entropy Production, and Stochastic Thermodynamics.

## Features

- **Onsager Relations**: Linear, Non-Linear, Stochastic, and Time-Dependent implementations
- **Entropy Production**: Local, rate calculations, time-dependent analysis, and optimization
- **Stochastic Thermodynamics**: Fluctuation theorems, Fokker-Planck equations, Langevin dynamics, and more
- **Minimal Dependencies**: Only requires `Plots` for visualization; all other functionality is manually implemented

## Installation

```julia
using Pkg
Pkg.add("NET")
```

## Quick Start

```julia
using NET

# Linear Onsager Relations
L = OnsagerMatrix([1.0 0.5; 0.5 2.0])
F = [1.0 2.0; 3.0 4.0]
J = compute_fluxes(L, F)

# Entropy Production
σ = compute_local_entropy_prozduction(J, F)
```

## Documentation

See the individual module documentation for detailed usage examples.

## License

MIT License - see [LICENSE](LICENSE) file for details.
