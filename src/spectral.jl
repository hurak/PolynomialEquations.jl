"""
    x = spectralfactor(a::Polynomial[,method=:roots])

Find a spectral factor `x` for a scalar univariate polynomial `a`.

For a scalar univariate polynomial `a` conjugate symmetric with respect to the imaginary axis, that is, `a=b̃b`, for some `b`, where `b̃(s)=b(-s)` (for complex polynomials the coefficients are complex-conjugated on top of substituting `-s` for `s`), find a polynomial `x` such that `x̃x=a` and `x` is Hurwitz stable (all roots in the left half plane).

Several methods are available: `:roots`, `:xxx`, `:xxx` (`xxx` is just a placeholder for methods to be added, the `:roots` methods is not a reliable method but is included as a "baseline" solution).

# Examples

```julia
x
```
"""
function spectralfactor(a::Polynomial)
    an = a[degree(a)]
    r = roots(a)
    m = map(x->isless(real(x),0.0),r)       # Currently ignoring roots on imaginary axis. TO DO!
    s = r[m]                                # Picking just the roots with negative real parts.
    x = fromroots(s,var=a.var)*sqrt(an)
    return x
end


"""
    x = spectralfactor(a::LaurentPolynomial[,method=:roots))

Find a spectral factor `x` for a scalar univariate Laurent polynomial `a`.

For a scalar univariate Laurent polynomial `a` conjugate symmetric with respect to the unit circle, that is, `a=b̃b`, for some `b`, where `b̃(z)=b(1/z)` (for complex polynomials the coefficients are complex-conjugated on top of substituting `1/z` for `z`), find a polynomial `x` such that `x̃x=a` and `x` is Schur stable (all roots in the unit disk).

Several methods are available: `:roots`, `:xxx`, `:xxx` (`xxx` is just a placeholder for methods to be added).

# Examples

```julia
x
```
"""
function spectralfactor(a::LaurentPolynomial)
    x=1
end
