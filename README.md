# PolynomialEquations

[![Build Status](https://github.com/hurak/PolynomialEquations.jl/workflows/CI/badge.svg)](https://github.com/hurak/PolynomialEquations.jl/actions)
[![Build Status](https://travis-ci.com/hurak/PolynomialEquations.jl.svg?branch=master)](https://travis-ci.com/hurak/PolynomialEquations.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/hurak/PolynomialEquations.jl?svg=true)](https://ci.appveyor.com/project/hurak/PolynomialEquations-jl)
[![Coverage](https://codecov.io/gh/hurak/PolynomialEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/hurak/PolynomialEquations.jl)
[![Coverage](https://coveralls.io/repos/github/hurak/PolynomialEquations.jl/badge.svg?branch=master)](https://coveralls.io/github/hurak/PolynomialEquations.jl?branch=master)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hurak.github.io/PolynomialEquations.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hurak.github.io/PolynomialEquations.jl/dev)

A Julia package for solving equations with univariate polynomials. Solvers for [linear (Diophantine) equations](https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine_equations) (also called [Bezout identity](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity)) in the ring of univariate polynomials and quadratic equations of special kind (also called spectral factorization) are (planned to be) provided.

It must be emphasized that the package does not offer functionality for solving nonlinear equations with the nonlinearities in the form of multivariate polynomials.

## Available solvers

The equations currently solved by the package are

```math
\begin{equation}
a(s)x(s)+b(s)y(s) = 0
\end{equation}
```
and

```math
\begin{equation}
a(s)x(s)+b(s)y(s) = c(s)
\end{equation}
```
where `a`, `b` and `c` are given univariate polynomials in the variable `s` and `x` and `y` are the polynomials to be found.

### Examples in Julia

```julia
julia> a = Polynomial([1,2,3],:s)
Polynomial(1 + 2*s + 3*s^2)

julia> b = Polynomial([4,5],:s)
Polynomial(4 + 5*s)

julia> c = Polynomial([1,1],:s)
Polynomial(1 + s)

julia> a = a*c
Polynomial(1 + 3*s + 5*s^2 + 3*s^3)

julia> b = b*c
Polynomial(4 + 9*s + 5*s^2)

julia> x, y = axby0(a,b)
Degree of x: 1  Dimension of nullspace of S: 1
Degree of x: 0  Dimension of nullspace of S: 0
(Polynomial(-0.5393598899705949 - 0.6741998624632413*s), Polynomial(0.1348399724926485 + 0.2696799449852973*s + 0.4045199174779448*s^2))

julia> a*x+b*y
Polynomial(-8.881784197001252e-16 - 4.440892098500626e-16*s - 8.881784197001252e-16*s^2 - 8.881784197001252e-16*s^3)
```

and

```julia
julia> a = Polynomial([1,2,3],:s)
Polynomial(1 + 2*s + 3*s^2)

julia> b = Polynomial([4,5],:s)
Polynomial(4 + 5*s)

julia> c = Polynomial([6,7,8],:s)
Polynomial(6 + 7*s + 8*s^2)

julia> x,y = axbyc(a,b,c)
(Polynomial(4.181818181818182), Polynomial(0.4545454545454547 - 0.909090909090909*s))

julia> a*x+b*y-c
Polynomial(8.881784197001252e-16*s)
```

## To do

The problems for which solvers are not yet implemented in this package are the symmetric and conjugate linear equations

```math
a(s)x(-s)+b(-s)y(s) = c(s)
```

```math
a(s)x(-s)+a(-s)x(s) = 2b(s)
```
```math
a(z)x(1/z)+b(1/z)y(z) = c(z)+d(1/z)
```
and

```math
a(z)x(1/z)+a(1/z)x(z) = b(z)+b(1/z)
```
in which `a`, `b` and `c` are given polynomials and `x` and `y` are the polynomials to be found.

A plan is to implement some solver(s) for the quadratic equations too, namely the equations of the type

```math
a(s)a(-s) = x(-s)x(s)
```
and

```math
a(z)a(1/z) = x(z)x(1/z)
```
where `a` is some given univariate polynomial and `x` is a Hurwitz stable and Schur stable polynomial, respectively, That is, all its roots are inside the open left half plane and an open unit circle, respectively.
