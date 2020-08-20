"""
    scale(a::Polynomial,ρ::Number)

Scale the univariate polynomial `a` with the variable `s` given by ``a(s) = a_0   + a_1 s + a_2 s^2 + ... + a_n s^n`` by the scalar positive real number `ρ`. This amounts to replacing the original variable `s` with a new scaled variable `ρs`.

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

Return the conjugate of a given univariate polynomial `a` given by ``a(s) = a_0 + a_1 s + a_2 s^2 + ... + a_n s^n`` with respect to the imaginary axis, that is, return the polynomial ``\\tilde a(s) = \\bar a(-s)= \\bar a_0 - \\bar a_1 s + \\bar a_2 s^2 + ... \\pm \\bar a_n s^n``.

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

Fof a polynomial ``a(s) = a_0 + a_1 s + a_2 s^2 + ... + a_n s^n``, return the polynomial ``r(s) = \\bar a_n + \\bar a_{n-1} s + \\bar a_{n-2} s^2 + ... + \\bar a_0 s^n``.

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

Return the conjugate polynomial for a given polynomial ``a(z) = a_0 + a_1 z + a_2 z^2 + ... + a_n z^n`` with respect to the unit circle in the complex plane, that is, return the Laurent polynomial ``\\tilde a(z) = \\bar a_n z^{-n} + \\bar a_{n-1} z^{-n+1} + \\bar a_{n-2} z^{-n+2} + ... + \\bar a_0``.

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
