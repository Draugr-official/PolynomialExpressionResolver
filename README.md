# Polynomial Expression Resolver

Polynomial Expression Resolver is a personal study of mathematics that helps you find a polynomial that fits an input, 2d graph.

# Features
- Polynomial interpolation by picking points and solving equations.
    - For low-degree polynomials, it builds a system of equations using a Vandermonde matrix and solves it using Gaussian elimination.
- If a fit isn't found, uses least-squares quadratic instead.
- In the last attempt, it’ll look for simple arithmetic formulas (like “2 * x + 3”).

# Notes
- It uses Gaussian elimination (with partial pivoting) for solving systems. in the case of a singular matrix, you’ll get null instead of hallucinations.
- Least-squares for quadratics uses a QR decomposition. It builds the right matrix and uses a modified Gram–Schmidt.
- Uses a tolerance (default 1e-7) to avoid issues with floating-point precision.

# Example
Input:
```
x-axis> 0, 1, 2, 3, 4, 5
y-axis> 30, 32, 34, 36, 38, 40
```

Result:
```
Polynomial fit found:
2*x + 30
```