"""
    x = spectralfactor(a::Polynomial[,method=:roots])

Find a spectral factor `x` for a scalar univariate polynomial `a`.

For a scalar univariate polynomial `a` conjugate symmetric with respect to the imaginary axis, that is, `a=b̃b`, for some `b`, where `b̃(s)=b(-s)` (for complex polynomials the coefficients are complex-conjugated on top of substituting `-s` for `s`), find a polynomial `x` such that `x̃x=a` and `x` is Hurwitz stable (all roots in the left half plane).

Several methods are available: `:roots`, `:xxx`, `:xxx` (`xxx` is just a placeholder for methods to be added, the `:roots` methods (based on computing the roots) is just a baseline solution not indended for serious use.

# Examples

```julia
julia> b = Polynomial([1.0,2.0,3.0],:s)*Polynomial([1,-1],:s)
Polynomial(1.0 + 1.0*s + 1.0*s^2 - 3.0*s^3)

julia> a = cconj(b)*b
Polynomial(1.0 + 1.0*s^2 + 7.0*s^4 - 9.0*s^6)

julia> roots(a)
6-element Array{Complex{Float64},1}:
  -1.0000000000000004 + 0.0im
 -0.33333333333333276 - 0.4714045207910314im
 -0.33333333333333276 + 0.4714045207910314im
  0.33333333333333326 - 0.47140452079103173im
  0.33333333333333326 + 0.47140452079103173im
   1.0000000000000024 + 0.0im

julia> x = spectralfactor(a)
Polynomial(0.9999999999999987 + 2.9999999999999956*s + 4.999999999999998*s^2 + 3.0*s^3)

julia> roots(x)
3-element Array{Complex{Float64},1}:
  -1.0000000000000022 - 2.220446049250313e-16im
 -0.33333333333333287 + 0.47140452079103157im
 -0.33333333333333254 - 0.4714045207910313im
```
"""
function spectralfactor(a::Polynomial)
    an = a[degree(a)]
    r = roots(a)
    m = map(x->isless(real(x),0.0),r)       # Currently ignoring roots on imaginary axis. TO DO!
    s = r[m]                                # Picking just the roots with negative real parts.
    x = fromroots(s,var=a.var)*sqrt(abs(an))
    return x
end


"""
    x = spectralfactor(a::LaurentPolynomial[,method=:roots))

Find a spectral factor `x` for a scalar univariate Laurent polynomial `a`.

For a scalar univariate Laurent polynomial `a` conjugate symmetric with respect to the unit circle, that is, `a=b̃b`, for some `b`, where `b̃(z)=b(1/z)` (for complex polynomials the coefficients are complex-conjugated on top of substituting `1/z` for `z`), find a polynomial `x` such that `x̃x=a` and `x` is Schur stable (all roots in the unit disk).

Several methods are available: `:roots`, `:xxx`, `:xxx` (`xxx` is just a placeholder for methods to be added). The method based on computing the roots is just a baseline solution not indended for serious use.

# Examples

```julia
julia> b = convert(LaurentPolynomial,fromroots([3,-1/2],var=:z));
julia> a = b*dconj(b)
LaurentPolynomial(-1.5*z⁻² + 1.25*z⁻¹ + 9.5 + 1.25*z - 1.5*z²)

julia> roots(a)
4-element Array{Float64,1}:
 -2.000000000000001
 -0.5
  0.33333333333333354
  2.999999999999998

julia> x = spectralfactor(a)
LaurentPolynomial(-0.5000000000000003 + 0.4999999999999994*z + 3.0*z²)

julia> roots(x)
2-element Array{Float64,1}:
 -0.5
  0.33333333333333354
```
"""
function spectralfactor(a::LaurentPolynomial)
    an = a[degree(a)]
    r = roots(a)
    m = map(x->isless(abs(x),1.0),r)    # Currently ignoring roots on the unit circle. TO DO!
    s = r[m]                            # Picking just the roots inside the unit disk.
    x = convert(LaurentPolynomial,fromroots(s,var=a.var))
    xx = x*dconj(x)                     # Finding α could be more efficent, should there be a need.
    α = sqrt(a[0]/(xx[0]))
    x = α*x
    return x
end
