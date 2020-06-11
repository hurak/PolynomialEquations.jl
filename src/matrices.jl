"""
    ltbtmatrix(a,c)

Build a lower triangular banded Toeplitz (LTBT) matrix corresponding to a univariate polynomial `a` and with the number of columns given by the integer `c`. The matrix is stored as dense.

# Examples

```julia
julia> a = Polynomial([1.0,2.0,3.0,4.0],:s)
julia> ltbtmatrix(a,3)
6×3 Array{Real,2}:
 1.0  0    0
 2.0  1.0  0
 3.0  2.0  1.0
 4.0  3.0  2.0
 0    4.0  3.0
 0    0    4.0
```
"""
function ltbtmatrix(a::Polynomial,c::Signed)
    da = degree(a)          # degree of the polynomial a
    r = da+c                # number of rows of the resulting matrix
    T = [i >= j && i <= (j+da) ? a[i-j] : 0.0 for i=1:r, j=1:c]
    return T
end

"""
    sylvestermatrix(a,b[,degx])

Build a Sylvester matrix for a pair of univariate polynomials `a` and `b`. The matrix is obtained by stacking horizontally two lower triangular banded Toeplitz matrices corresponding to the two polynomials. If the third parameter `degx` is not specified, it is set to `degree(b)`. The matrix is stored as dense.

# Examples

```julia
julia> a = Polynomial([1,2,3],:s)
Polynomial(1 + 2*s + 3*s^2)

julia> b = Polynomial([4,5],:s)
Polynomial(4 + 5*s)

julia> S = sylvestermatrix(a,b)
3×3 Array{Float64,2}:
 1.0  4.0  0.0
 2.0  5.0  4.0
 3.0  0.0  5.0
 ```
"""
function sylvestermatrix(a::Polynomial,b::Polynomial;degx=degree(b)-1)
    da = degree(a)
    db = degree(b)
    dx = degx
    dy = da+dx-db
    n = 1+da+dx                    # number of rows (default: da+db)
    m = 1+dx+1+dy                  # number of columns (default: da+db)
    S = zeros(Float64,n,m)         # "hard-wired" to Float64 because we cannot exploit the knowledge of the data type anyway. TO DO: should be type-parameterized.
    for i = 1:n
        for j = 1:m
            if j <= dx+1           # if in the "left part" of the Sylvester matrix
                S[i,j] = i >= j && i <= (j+da) ? a[i - j] : 0.0
            else
                k = j-(dx+1)       # index of the column in the "right part" of the Sylvester matrix
                S[i,j] = i >= k && i <= (k+db) ? b[i - k] : 0.0
            end
        end
    end
    return S
end
