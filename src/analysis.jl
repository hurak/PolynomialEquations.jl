"""
    iscoprime(a,b)

Check if the pair of polynomials `a` and `b` is coprime, that is, it has no common factors.

## Examples

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
    S = sylvestermatrix(a,b,degx=dx)
    return rank(S)==(da+db)
end
