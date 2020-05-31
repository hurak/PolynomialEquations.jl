"""
    x,y = axby0(a,b)

Solve the homogeneous linear Diophantine equation ``ax+by=0`` with scalar univariate polynomials.

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
    #println("Degree of x: $dx  Dimension of nullspace of S: $ds")
    while !isempty(xy)
        x = Polynomial(xy[1:dx+1],a.var)
        y = Polynomial(xy[dx+2:end],a.var)
        dx = dx-1         # simply reducing the anticipated degrees of the solutions
        S = sylvestermatrix(a,b,degx=dx)
        xy = nullspace(S)
        ds = size(xy,2)
        #println("Degree of x: $dx  Dimension of nullspace of S: $ds")
    end
    return (x,y)
end

"""
   x,y = axbyc(a,b,c)

Solve the linear Diophantine equation ``ax+by=c`` with univariate polynomials.

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

"""
    axycminl1(a,c[,dymax=100,emax=1e-4])

Solve the linear Diophantine equation ``ax+y=c`` for given univariate polynomials `a` and `c`, that is, find two univariate polynomials `x` and `y` satisfying the equation and, moreover, with the 1-norm of the vector of coefficients of the polynomial `y` minimized.

The degree of the minimizing polynomial `y` is not known apriori. It is, however, known that such minimizing polynomials of a finite degree exists. The solver then starts with smalest possible degree of `y` and goes up to `dmax` while checking if the improvement with respect to the previous iteration is better than `emax`.

# Examples

```julia
julia> a = Poly([-0.12, -0.29, 0.71, 2.79, 2.92, 1],:d)
Poly(-0.12 - 0.29*d + 0.71*d^2 + 2.79*d^3 + 2.92*d^4 + 1.0*d^5)

julia> c = Poly([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:d)
Poly(1.86 - 0.34*d - 1.14*d^2 - 0.21*d^3 + 1.19*d^4 - 1.12*d^5)

julia> x, y, val = axycminl1(a,c);
```
Verify that the solutions are truly satisfying the original equation
```julia
julia> a*x+y≈c
```
Note that the degree of the resulting polynomial `y` might look high but it is only an artefact of the numerical solution.
```julia
julia> degree(chop(y,atol=1e-3))
```

# References

Hurak, Z., A. Böttcher, and M. Sebek, "Minimum Distance to the Range of a Banded Lower Triangular Toeplitz Operator in l1 and Application in l1-Optimal Control", SIAM Journal on Control and Optimization, vol. 45, issue 1: SIAM, pp. 107-122, 2006. [DOI: 10.1137/S0363012903437940](http://dx.doi.org/10.1137/S0363012903437940).
"""
function axycminl1(a::Polynomial,c::Polynomial;dymax=100,emax=1e-4)
    da = degree(a)                      # degree of a
    dc = degree(c)                      # degree of c
    dymin = max(da,dc)                  # initial degree of y
    b = append!(copy(coeffs(c)),zeros(dymin-dc)) # initial length of the right hand column
    l1y_pre = Inf                       # initial optimal value (l1 norm of the coefficients of y)
    for dy = dymin:dymax
        dx = dy-da
        T = ltbtmatrix(a,dx+1)          # lower triangular banded Toeplitz matrix
        model = Model(Clp.Optimizer)
        @variable(model, x[1:dx+1])     # coefficients of the x polynomial
        @variable(model, y⁺[1:dy+1]>=0) # coefficients of the y polynomial, the plus part
        @variable(model, y⁻[1:dy+1]>=0) # coefficients of the y polynomial, the minus part
        @constraint(model, con, T*x + y⁺ - y⁻ .== b)
        @objective(model, Min, sum(y⁺+y⁻))
        set_silent(model::Model)
        optimize!(model)
        l1y = objective_value(model)
        xₛ = value.(x)
        y⁺ₛ = value.(y⁺)
        y⁻ₛ = value.(y⁻)
        yₛ = y⁺ₛ-y⁻ₛ
        @printf("|y|₁=%f \t\t Deg(x)=%d \t Deg(y)=%d\n",l1y,dx,dy)   # print the value just for diagnostics
        b = append!(b,0)                # extend the b for the next iteration
        global xₛ
        global yₛ
        global l1y
        e = l1y_pre-l1y                 # diff. between the value of the current and the previous
        abs(e)<emax && break            # if not improved enough, break
        l1y_pre = l1y
    end
    x = Polynomial(xₛ,a.var)
    y = Polynomial(yₛ,a.var)
    return (x,y,l1y)
end
