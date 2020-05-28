"""
    x,y = axby0(a,b)

Solve the homogeneous linear (Diophantine) equation ``ax+by=0`` with scalar univariate polynomials.

The problem corresponds to finding the minimum-order representation ``x/y`` of the fraction ``b/a``, that is, it solves ``b/a = x/y`` for a given fraction ``b/a`` such that `y` is of minimum possible degree.

# Examples

```julia
julia> a = Polynomial([1,2,3],:s); b = Polynomial([4,5],:s); c = Polynomial([1,1],:s);
julia> a = a*c; b = b*c
julia> x, y = axby0(a,b)
(Polynomial(-0.5393598899705949 - 0.6741998624632413*s), Polynomial(0.1348399724926485 + 0.2696799449852973*s +
0.4045199174779448*s^2))
```
"""
function axby0(a::Polynomial,b::Polynomial)
    da = degree(a)
    db = degree(b)
    rd = da-db            # relative degree of the transfer function b/a
    dx = db-1             # no need to start with deg(x) = deg(b) because such solution is guaranteed to exist:
    x = -b
    y = a
    S = sylvestermatrix(a,b,degx=dx)
    xy = nullspace(S)
    ds = size(xy,2)
    println("Degree of x: $dx  Dimension of nullspace of S: $ds")
    while !isempty(xy)
        x = Polynomial(xy[1:dx+1],a.var)
        y = Polynomial(xy[dx+2:end],a.var)
        dx = dx-1         # simply reducing the anticipated degrees of the solutions
        S = sylvestermatrix(a,b,degx=dx)
        xy = nullspace(S)
        ds = size(xy,2)
        println("Degree of x: $dx  Dimension of nullspace of S: $ds")
    end
    return (x,y)
end

"""
   x,y = axbyc(a,b,c)

Solve the linear (Diophantine) equation ``ax+by=c`` with univariate polynomials.

# Examples

```jldoctest
julia> a = Polynomial([1,2,3],:s);  b = Polynomial([5,6],:s); c = Polynomial([6,7,8],:s);

julia> x, y = axbyc(a,b,c)
(Polynomial(4.181818181818182), Polynomial(0.4545454545454546 - 0.9090909090909091*s))
```
"""
function axbyc(a::Polynomial,b::Polynomial,c::Polynomial)
    da = degree(a)
    db = degree(b)
    dc = degree(c)
    dx = db - 1
    dy = da -1
    S = sylvestermatrix(a,b)
    xy = S\coeffs(c)
    x = xy[1:dx+1]
    y = xy[dx+2:(da+db)]
    x = Polynomial(x,a.var)
    y = Polynomial(y,a.var)
    return (x,y)
end
