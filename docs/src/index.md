```@meta
CurrentModule = PolynomialEquations
```

# PolynomialEquations

A Julia package for solving equations with univariate polynomials. Solvers for [linear (Diophantine) equations](https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine_equations) (also called [Bezout identity](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity)) in the ring of univariate polynomials and quadratic equations of special kind (also called spectral factorization) are (planned to be) provided.

!!! note "Equations with polynomials and not polynomial equations in real variables"
    The package does not offer functionality for solving nonlinear equations (with real variables), in which the nonlinearities come in the form of multivariate polynomials (in those real variables). Such functionality can be found elsewhere.

## Linear (Diophantine) equations

```@docs
axbyc
axby0
axb
axaxbb
axbycd
axycminl1
```

## Quadratic equations — spectral factorization

Nothing yet.

## Utilities

### Useful matrices built from the coefficients of polynomials

```@docs
ltbtmatrix
sylvestermatrix
```
### Some operations on polynomials

```@docs
scale
cconj
dconj
conjreciprocal
```
### Analysis of polynomials

```@docs
iscoprime
ishurwitzstable
isschurstable
```
