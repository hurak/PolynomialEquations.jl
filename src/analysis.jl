"""
    iscoprime(a,b)

Check if the pair of polynomials `a` and `b` is coprime, that is, it has no common factors.

# Examples

```julia
julia> a = Polynomial([1,2,3],:s);
julia> b = Polynomial([1,2],:s);
julia> iscoprime(a,b)
true

julia> a = Polynomial([1,2,3],:s);
julia> b = a*Polynomial([1,2],:s);
julia> iscoprime(a,b)
false
```
"""
function iscoprime(a::Polynomial,b::Polynomial)
    da = degree(a)
    db = degree(b)
    dx = db-1
    S = sylvesterresultantmatrix(a,b)
    return rank(S)==(da+db)
end

"""
    isschurstable(a)

Check if the polynomial `a` is Schur stable, that is, it has all its roots inside the unit disk.

# Examples

```julia
julia> a = fromroots([-1/2,0, 2/3,-1/3im, 1/3im]);
julia> isschurstable(a)
true

julia> a = fromroots([2,0, 2/3,-3im, 3im]);
julia> !isschurstable(a)
false
```
"""
function isschurstable(a::Polynomial)
    r = roots(a)
    all(abs.(r) .<= 1)
end

"""
    ishurwitzstable(a::Polynomial)

Check if the `a` is Schur stable, that is, it has all its roots inside the left complex half-plane.

# Examples

```julia
julia> a = fromroots([-1, -2, -3+4im, -3-4im]);
julia> ishurwitzstable(a)
true

julia> a = fromroots([-1, -2, 3+4im, 3-4im]);
julia> !ishurwitzstable(a)
false
```
"""
function ishurwitzstable(a::Polynomial)
    r = roots(a)
    return all(real(r) .<= 0)
end

"""
    isconjsymmetric(a::Polynomial)

Check if the scalar univariate polynomial `a` is conjugate symmetric, that is, if `a=ã`.

For a scalar univariate polynomial `a` with real coefficients, check if `a(s)=a(-s)`. If the coefficients are complex, the condition then extends to `a(s)=ā(-s)`, where `ā` is obtained from `a` by complex-conjugating the coefficients. A necessary and sufficient condition for conjugate symmetry is that the coefficients with odd powers of the variable are zero and the coefficients with even powers are real.

# Examples
```julia
julia> a = Polynomial([1,0,2,0,3],:s)
Polynomial(1 + 2*s^2 + 3*s^4)

julia> isconjsymmetric(a)
true
```
"""
function isconjsymmetric(a::Polynomial)
    return isequal(a,cconj(a))
end

"""
    isconjsymmetric(a::LaurentPolynomial)

Check if the scalar univariate Laurent polynomial `a` is conjugate symmetric, that is, if `a=ã`.

For a scalar univariate Laurent polynomial `a` with real coefficients, check if `a(z)=a(1/z)`. If the coefficients are complex, the condition then extends to `a(z)=ā(1/z)`, where `ā` is obtained from `a` by complex-conjugating the coefficients. A necessary and sufficient condition is that the coefficients are real and their order in the coefficient vector can be flipped with no impact.

# Examples
```julia
julia> a = LaurentPolynomial([3,2,1,2,3],-2:2,:z)
LaurentPolynomial(3*z⁻² + 2*z⁻¹ + 1 + 2*z + 3*z²)

julia> isconjsymmetric(a)
true
```
"""
function isconjsymmetric(a::LaurentPolynomial)
    return isequal(a,dconj(a))
end
