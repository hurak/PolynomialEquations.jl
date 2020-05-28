"""
    scale(a,ρ)

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
