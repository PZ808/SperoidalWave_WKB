# Gauss–Manin / Picard–Fuchs computation in Julia

This project implements a symbolic Julia/OSCAR pipeline for computing Gauss–Manin connection matrices and scalar Picard–Fuchs equations for the genus-one family
$$
f(x,y,z;E,\mu_0,\mu_1,\mu_2)
=
\frac12 y^2
+x^4
-\mu_2 x^3 z
+(E-2)x^2 z^2
+(\mu_1+\mu_2)x z^3
+(\mu_0-E+1)z^4
$$
in weighted projective space $\mathbb P^{(1,2,1)}$.

The main goal is to compute the variation of period integrals with respect to the parameters
$$
(E,\mu_0,\mu_1,\mu_2),
$$
using Griffiths–Dwork reduction, Gröbner basis methods, and the resulting Gauss–Manin system.

---

## Mathematical setup

We work with the two-dimensional de Rham basis
\[
\omega_1=\frac{\Omega}{f},
\qquad
\omega_2=\frac{z^4\Omega}{f^2},
\]
where \(\Omega\) is the standard weighted projective ambient form.

For each parameter \(\lambda \in \{E,\mu_0,\mu_1,\mu_2\}\), the Gauss–Manin matrix is defined by
\[
\partial_\lambda \Pi = A_\lambda \Pi,
\]
where \(\Pi\) is the period vector and \(A_\lambda\) is a \(2\times 2\) matrix of rational functions in the parameters.

From the \(E\)-sector matrix \(A_E\), one can eliminate one component and recover scalar second-order Picard–Fuchs equations.

---

## What this code does

The code provides tools to:

- build the weighted polynomial ring and Jacobian ideal
- compute Gröbner bases and transformation matrices
- compute lifts with respect to the Jacobian ideal efficiently
- perform Griffiths–Dwork pole reduction
- project reduced classes onto the basis \(\{\Omega/f,\ z^4\Omega/f^2\}\)
- construct the Gauss–Manin matrices \(A_E, A_{\mu_0}, A_{\mu_1}, A_{\mu_2}\)
- check flatness / Picard–Fuchs consistency
  $$
  \partial_{\lambda_1}A_{\lambda_2}
    -
  \partial_{\lambda_2}A_{\lambda_1}
    +
  [A_{\lambda_1},A_{\lambda_2}]
  =0
  $$
- extract scalar Picard–Fuchs equations from the \(E\)-sector

---

## Main implementation idea

A major bottleneck in a naive Griffiths–Dwork implementation is the repeated computation of **lifts**
$$
p=\sum_i a_i J_i
$$
with respect to the Jacobian generators \(J_i\).

Instead of solving this linear problem from scratch each time, this code uses OSCAR's

- `groebner_basis_with_transformation_matrix`
- `reduce_with_quotients`

to precompute a Gröbner basis together with the transformation matrix back to the original generators. 
This makes the lift step substantially faster and was essential for making the full Gauss–Manin computation practical.

---

## Requirements

- Julia 1.10+
- [OSCAR.jl](https://oscar-system.github.io/Oscar.jl/stable/)

Depending on your setup, you may also want IJulia or a Jupyter notebook environment.

---

## Installation

Clone the repository and instantiate the Julia environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()