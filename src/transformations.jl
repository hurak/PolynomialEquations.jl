"""
    scale(a::Polynomial,ρ::Number)

Scale the univariate polynomial `a` with the variable `s` given by ``a(s) = a₀ + a₁ s + a₂ s² + … + aₙ sⁿ`` by the scalar positive real number `ρ`. This amounts to replacing the original variable `s` with a new scaled variable `ρs`.

# Examples

```julia
julia> a = Polynomial(1:5,:s)
julia> ρ=10
julia> scale(a,ρ)
Polynomial(1 + 20*s + 300*s^2 + 4000*s^3 + 50000*s^4)
```
"""
function scale(a::Polynomial,ρ::Number)
    c = deepcopy(coeffs(a))     # coefficent vector of the original polynomial
    for k=1:degree(a)           # k is the power to which s is raised in s^k
        c[k+1] *= ρ^k           # coefficient vector of the new scaled polynomial
    end
    p = Polynomial(c,a.var)     # converting the coefficient vector to a polynomial
end

"""
    cconj(a::Polynomial)

Return the conjugate of a given univariate polynomial `a` given by ``a(s) = a₀ + a₁ s + a₂ s² + … + aₙ sⁿ`` with respect to the imaginary axis, that is, return the polynomial ``̃ã(s) = ̄ā(-s)=ā₀ - ̄ā₁s + ā₂ s² + … ± ̄āₙ sⁿ``.

It is used in the analysis and synthesis of continuous-time filters and controllers. This is reflected in prepending the letter `c` to the `conj` function name.

# Examples

```julia
julia> a = Polynomial(1:5,:s)
Polynomial(1 + 2*s + 3*s^2 + 4*s^3 + 5*s^4)

julia> p = cconj(a)
Polynomial(1 - 2*s + 3*s^2 - 4*s^3 + 5*s^4)

julia> a = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
Polynomial((1 + 1im) + (2 + 2im)*s + (3 + 3im)*s^2 + (4 + 4im)*s^3 + (5 + 5im)*s^4)

julia> p = cconj(a)
Polynomial((1 - 1im) - (2 - 2im)*s + (3 - 3im)*s^2 - (4 - 4im)*s^3 + (5 - 5im)*s^4)
```
"""
function cconj(a::Polynomial)
    c = conj(a)
    cc = deepcopy(coeffs(c))
    cc[2:2:end] = -cc[2:2:end]
    return Polynomial(cc,a.var)
end

"""
    conjreciprocal(a::Polynomial)

Return the conjugate reciprocal polynomial for a given polynomial `a`.

Fof a polynomial ``a(s) = a₀ + a₁ s + a₂ s² + … + aₙ sⁿ``, return the polynomial ``r(s) = ̄āₙ + ̄āₙ₋₁ s + ̄āₙ₋₂ s² + … + ̄ā₀ sⁿ``.

# Examples

```julia
julia> a = Polynomial(1:5,:s)
Polynomial(1 + 2*s + 3*s^2 + 4*s^3 + 5*s^4)

julia> p = conjreciprocal(a)
Polynomial(5 + 4*s + 3*s^2 + 2*s^3 + s^4)

julia> c = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
Polynomial((1 + 1im) + (2 + 2im)*s + (3 + 3im)*s^2 + (4 + 4im)*s^3 + (5 + 5im)*s^4)

julia> p = conjreciprocal(c)
Polynomial((5 - 5im) + (4 - 4im)*s + (3 - 3im)*s^2 + (2 - 2im)*s^3 + (1 - 1im)*s^4)
```
"""
function conjreciprocal(a::Polynomial)
    c = conj(a)
    cc = deepcopy(coeffs(c))
    pc = cc[end:-1:1]
    p = Polynomial(pc,a.var)
end

"""
    dconj(a::LaurentPolynomial)

Return the conjugate polynomial for a given polynomial ``a(z) = a₀ + a₁ z + a₂ z² + ... + aₙ zⁿ`` with respect to the unit circle in the complex plane, that is, return the Laurent polynomial ``ã(z) = āₙ z⁻ⁿ + āₙ₋₁ z⁻ⁿ⁺¹ + āₙ₋₂ z⁻ⁿ⁺² + … + ā₀``.

The function is only defined for `LaurentPolynomial` type even if it is used to represent a standard polynomial (with no negative powers).

It is used in the analysis and synthesis of discrete-time filters and controllers. This is reflected in prepending the letter `d` to the `conj` function name.

# Examples

```julia
julia> a = LaurentPolynomial(1:5,:z)
LaurentPolynomial(1 + 2*z + 3*z² + 4*z³ + 5*z⁴)

julia> p = dconj(a)
LaurentPolynomial(5*z⁻⁴ + 4*z⁻³ + 3*z⁻² + 2*z⁻¹ + 1)

julia> c = LaurentPolynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:z)
LaurentPolynomial((1 + 1im) + (2 + 2im)*z + (3 + 3im)*z² + (4 + 4im)*z³ + (5 + 5im)*z⁴)

julia> p = dconj(c)
LaurentPolynomial((5 - 5im)*z⁻⁴ + (4 - 4im)*z⁻³ + (3 - 3im)*z⁻² + (2 - 2im)*z⁻¹ + (1 - 1im))
```
"""
function dconj(a::LaurentPolynomial)
    ac = deepcopy(coeffs(a))
    dc = conj.(ac[end:-1:1])             # flipping and conjugating the coeffs
    fr = -a.n[]:-a.m[]                   # flipping and negating the range
    d = LaurentPolynomial(dc,fr,a.var)
end

"""
    shift(p::LaurentPolynomial[,k::Integer=1])

Increase the powers of the variable in the Laurent polynomial by 1 or a given (possibly negative) number.

For a univariate Laurent polynomial `p(z) = pₘ zᵐ + p₋₁ z⁻¹ +  p₀ + p₁ z + … + pₙ zⁿ` and a given integer `k`, return a Laurent polynomial  `p(z) = pₘ zᵐ⁺ᵏ + p₋₁ z⁻¹⁺ᵏ +  p₀ zᵏ + p₁ zᵏ⁺¹ + … + pₙ zⁿ⁺ᵏ`. If `k` is not specified, it is assumed that `k=1`.

# Examples

```julia
julia> p = LaurentPolynomial([1,2,3,4],-2:1,:z)
LaurentPolynomial(z⁻² + 2*z⁻¹ + 3 + 4*z)

julia> shift(p)
LaurentPolynomial(z⁻¹ + 2 + 3*z + 4*z²)

julia> shift(p,-3)
LaurentPolynomial(z⁻⁵ + 2*z⁻⁴ + 3*z⁻³ + 4*z⁻²)
```
"""
function shift(p::LaurentPolynomial,k::Integer=1)
    r = range(p)
    s = LaurentPolynomial(coeffs(p),r.+k,p.var)
    return s
end
