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

### Homogeneous linear diophantine equation
```math
a(s)x(s)+b(s)y(s) = 0
```
where `a` and `b` are given univariate polynomials in the variable `s` and `x` and `y` are the polynomials to be found.

#### Example
```julia
julia> a = Polynomial([1,2,3],:s);
julia> b = Polynomial([4,5],:s);
julia> c = Polynomial([1,1],:s);
julia> a = a*c;
julia> b = b*c;

julia> x, y = axby0(a,b)
(Polynomial(-0.5393598899705949 - 0.6741998624632413*s), Polynomial(0.1348399724926485 + 0.2696799449852973*s + 0.4045199174779448*s^2))

julia> a*x+b*y
Polynomial(-8.881784197001252e-16 - 4.440892098500626e-16*s - 8.881784197001252e-16*s^2 - 8.881784197001252e-16*s^3)
```

### Linear diophantine equation
```math
a(s)x(s)+b(s)y(s) = c(s)
```
where `a`, `b` and `c` are given univariate polynomials in the variable `s` and `x` and `y` are the polynomials to be found.

#### Example
```julia
julia> a = Polynomial([1,2,3],:s);
julia> b = Polynomial([4,5],:s);
julia> c = Polynomial([6,7,8],:s);

julia> x,y = axbyc(a,b,c)
(Polynomial(4.181818181818182), Polynomial(0.4545454545454547 - 0.909090909090909*s))

julia> a*x+b*y≈c
true
```

### Symmetrix conjugate equation (continuous-time case)  

```math
a(s)x(-s)+a(-s)x(s) = b(s)+b(-s)
```

where `a`, `b` and `c` are given univariate polynomials in the variable `s` and `x` and `y` are the polynomials to be found.

If the coefficients of the polynomials are complex, then not only `-s` is substituted for `s` but also the coefficients are complex-conjugated.

Formally, both real and complex cases are captured in the following notation

```math
a(s)x̃(s)+ã(s)x(s) = b(s)+b̃(s)
```
where tilde denotes *conjugation of a polynomial with respect to the imaginary axis*, that is, the conjugate polynomial evaluates to a complex conjugate value of the original polynomial on the imaginary axis.  

#### Example
```julia
julia> a = Polynomial([-0.12, -0.29, 1],:s);
julia> b = Polynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:s);

julia> x = axaxbb(a,b)
Polynomial(-15.50000000000003 + 50.0096551724139*s + 1.19*s^2)

julia> cconj(a)*x+a*cconj(x)≈b+cconj(b)
true
```

### Nonsymmetric conjugate equation (continuous-time case)
As a generalization of the previous symmetric conjugate equation, there is also a solver for a nonsymmetric conjugation.

```math
a(s)x(-s)+b(-s)y(s) = c(s)+d(-s)
```
#### Example
```julia
julia> a = Polynomial([-0.12, -0.29, 1],:s);
julia> b = Polynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:s);
julia> c = Polynomial([1.2, 3.4, 5.6],:s);
julia> d = Polynomial([3.4, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2],:s);

julia> x,y = axbycd(a,b,c,d)
(Polynomial(13.57846778138169 + 9.81569937794185*s + 1.6765177065455532*s^2 + 12.031635021058513*s^3 + 1.5485480613952116*s^4), Polynomial(3.349148459013872 - 0.31120362624572473*s))

julia> a*cconj(x)+cconj(b)*y≈c+cconj(d)
true
```

## To do

The symmetric and nonsymmetric conjugate equations for Laurent polynomials

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

## Related Julia packages
- [Polynomials](https://github.com/JuliaMath/Polynomials.jl) – provides the basic data types – `Polynomial` and `LaurentPolynomial`.
- [PolynomialMatrices](https://github.com/JuliaPolynomialMatrices/PolynomialMatrices.jl) – a package defining a new data type – a polynomial matrix (or matrix polynomial). In principle, many if not all equations with polynomials can be reformulated as equations with polynomial matrices. For example, the equation ``ax+by=0`` can be formulated as ``[a b][x;y]=0``, that is a polynomial nullspace problem ``Ax=0``, or similarly ``ax+by=c`` can be rewritten as ``AX =c``, therefore developing dedicated scalar versions of solvers might be viewed as redundant. Still, in many projects scalar polynomials are used exclusively (SISO filter and controller design) and avoiding the complications with polynomial matrices is useful.   
- [ControlSystems](https://github.com/JuliaControl/ControlSystems.jl) – a package for control system design. This package and the next  might be finally benefitting from the `PolynomialEquations` package. Some solver for linear diophantine equation with univariate polynomials is already implemented in `ControlSystems` package, namely the [dab](http://juliacontrol.github.io/ControlSystems.jl/latest/lib/synthesis/#ControlSystems.dab) function, but others are missing.
- [DSP](https://github.com/JuliaDSP/DSP.jl) – a package for digital signal processing.

## Related other software packages

- [Polynomial Toolbox for Matlab](http://polyx.com/)  - one of the most complete packages for numerical computations with univariate polynomials. Commercial.
