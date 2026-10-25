# Assignment 3 - SPH Kernel Documentation

## Overview
This document summarizes the implementation of the cubic spline kernel `W(x_i - x_j; h)` and its gradient `∇W` for the learnSPH framework.

## Kernel Function

### Definition
Cubic spline kernel (Monaghan):
```
W(x - x'; h) = (1/h^3) * s(q),  q = |x - x'| / h
s(q) = alpha * {
  2/3 - q^2 + 0.5*q^3,        0 <= q < 1
  (1/6)*(2-q)^3,              1 <= q < 2
  0,                          q >= 2
}
alpha = 3/(2*pi)
```

### Implementation
Files:
- `framework/learnSPH/kernel.h`: declarations and inline helpers
- `framework/learnSPH/kernel.cpp`: implementations of `cubic_spline(q)` and `cubic_grad_spline(q)`

Key points:
- `cubic_spline(q)` enforces `q >= 0` via assert. Returns `alpha * s(q)` form from the formula scaled correctly for all branches (branch 2 includes `alpha`).
- `cubic_spline(xi, xj, h)` computes `q = |xi-xj|/h` and returns `cubic_spline(q) / (h^3)`.
- Compact support: returns 0 for `q >= 2`.

## Kernel Gradient

### Analytical Form
For `W = (1/h^3) s(q)`, `q = r/h`, `r = |x_i - x_j|`:
- `∇_{x_i} q = (x_i - x_j)/(h r) = r/(h^2 q)` (with care at r=0)
- `∇_{x_i} W = (1/h^3) s'(q) ∇_{x_i} q = s'(q)/(h^4 dist) * (x_i - x_j)`

Derivative `s'(q)`:
- `0 <= q < 1`: `s'(q) = alpha * (-2q + 1.5*q^2)`
- `1 <= q < 2`: `s'(q) = -alpha/2 * (2-q)^2`
- `q >= 2`: `0`

Implementation: `cubic_grad_spline(q)` computes `s'(q)`, `cubic_spline_grad(xi,xj,h)` returns `s'(q)/(h^4 * dist) * r` for `dist > 1e-18` and `q < 2`, else zero vector.

### Special Case r = 0
When `x_i == x_j`, distance is 0. The gradient is returned as zero vector to avoid division by zero. Analytically the cubic spline kernel has a smooth maximum at `q=0` and is radial; the gradient of a radial function `W(r)` at `r=0` is 0 (by symmetry / odd derivative). So returning `Eigen::Vector3d::Zero()` is correct.

## Testing

Tests in `framework/tests/kernel_tests.cpp` (built by `framework/tests/CMakeLists.txt`).

### Branch Coverage (a)
- Branch 1 `[0,1)`: `q=0` -> `alpha*2/3`, `q=0.5` -> formula check, continuity from below at `q->1^-`
- Branch 2 `[1,2)`: `q=1.5` -> `alpha/6*(2-q)^3`, continuity at `q=1` between branches
- Branch 3 `[2,∞)`: `q>=2` returns 0, tested with fixed values and random samples in `[2,100]`
- `W(xi,xj,h)` formula validated via `cubic_spline(q)/h^3`

### Gradient Tests (b),(c)
- Finite difference: `∇_x W ≈ (1/(2ε))(W(x+εe_i)-W(x-εe_i))` with `ε=1e-6`, checked for branch 1 and branch 2 configurations. All components match within `1e-6`.
- `r==0`: gradient is zero vector.
- Compact support: `q>2` gives zero gradient.

### Running Tests
```bash
cmake --build /home/cube/projects/SPH/SPH_Sim/build --target learnSPH_tests -j
./home/cube/projects/SPH/SPH_Sim/build/tests/learnSPH_tests "[kernel]"
./home/cube/projects/SPH/SPH_Sim/build/tests/learnSPH_tests "[grad]"
```

## Design Notes (d) - Performance Considerations
For constant smoothing length `h` in a simulation:
- Precompute constants once: `inv_h3 = 1/h^3`, `inv_h4 = 1/h^4`, `h2=h^2`, `h4=h^4` (avoid repeated divisions/powers)
- Avoid `sqrt` when possible: compare `r^2 < 4*h^2` for support test, compute `r` only when inside support
- Prefer multiplications over `pow` for small integer powers (e.g. `(2-q)^2`, `(2-q)^3`)
- Minimize memory indirection; pass vectors by const ref
- Branch prediction friendly checks (`q < 1`, `q < 2`)
- Keep simple API for testing; optimized version can take precomputed values

Current implementation prioritizes clarity/correctness for testing; optimized variants can be added when integrating into the simulation loop.
