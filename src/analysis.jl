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

"""
    isschurstable(a)

Check if the `a` is Schur stable, that is, it has all its roots inside the unit disk.

## Examples

```julia
julia> a = fromroots([-1/2,0, 2/3,-1/3im, 1/3im])
julia> isschurstable(a)
true

julia> a = fromroots([2,0, 2/3,-3im, 3im])
julia> !isschurstable(a)
false
```
"""

function isschurstable(a::Polynomial)
    r = roots(a)
    all(abs.(r) .<= 1)
end

"""
    ishurwitzstable(a)

Check if the `a` is Schur stable, that is, it has all its roots inside the left complex half-plane.

## Examples

```julia
julia> a = fromroots([-1, -2, -3+4im, -3-4im])
julia> ishurwitzstable(a)
true

julia> a = fromroots([-1, -2, 3+4im, 3-4im])
julia> !ishurwitzstable(a)
false
```
"""

function ishurwitzstable(a::Polynomial)
    r = roots(a)
    return all(real(r) .<= 0)
end
