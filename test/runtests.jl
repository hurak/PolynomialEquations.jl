using Polynomials
using PolynomialEquations
using Test
using LinearAlgebra

@testset "Scaling of a polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        ρ=10
        as = scale(a,ρ)
        r = Polynomial([1,20,300,4000,50000],:s)
        isequal(as,r)
    end
end

@testset "Conjugate of a polynomial with respect to the imaginary axis" begin
    @test begin
        a = Polynomial(1:5,:s)
        p = cconj(a)
        r = Polynomial([1,-2,3,-4,5],:s)
        isequal(p,r)
    end
    @test begin
        a = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
        p = cconj(a)
        r = Polynomial([1-1im, -(2-2im), 3-3im, -(4-4im), 5-5im],:s)
        isequal(p,r)
    end
end

@testset "Conjugate reciprocal polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        p = conjreciprocal(a)
        r = Polynomial([5,4,3,2,1],:s)
        isequal(p,r)
    end
    @test begin
        c = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
        p = conjreciprocal(c)
        r = Polynomial([5-5im, 4-4im, 3-3im, 2-2im, 1-1im],:s)
        isequal(p,r)
    end
end

@testset "Conjugate of a polynomial with respect to the unit circle" begin
    @test begin
        a = LaurentPolynomial(1:5,:z)
        d = dconj(a)
        r = LaurentPolynomial(5:-1:1,-4:0,:z)
        isequal(d,r)
    end
    @test begin
        a = LaurentPolynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:z)
        d = dconj(a)
        r = LaurentPolynomial([5-5im, 4-4im, 3-3im, 2-2im, 1-1im],-4:0,:z)
        isequal(d,r)
    end
end

@testset "Building a lower triangular banded Toeplitz matrix for a polynomial" begin
    @test begin
        a = Polynomial([1.0,2.0,3.0,4.0],:s)
        A = ltbtmatrix(a,3)
        R = [1.0  0    0;
             2.0  1.0  0;
             3.0  2.0  1.0;
             4.0  3.0  2.0;
             0    4.0  3.0;
             0    0    4.0]
        isequal(A,R)
    end
    @test begin
        a = Polynomial([1,2,3,4],:s)
        A = ltbtmatrix(a,3)
        R = [1  0  0;
             2  1  0;
             3  2  1;
             4  3  2;
             0  4  3;
             0  0  4]
        isequal(A,R)
    end
end

@testset "Building a lower triangular banded Toeplitz matrix for a Laurent polynomial" begin
    @test begin
        a = LaurentPolynomial([1,2,3,4,5],-2:2,:z)
        A = ltbtmatrix(a,3)
        R = [1  0  0;
             2  1  0;
             3  2  1;
             4  3  2;
             5  4  3;
             0  5  4;
             0  0  5]
        isequal(A,R)
    end
end

@testset "Building a Sylvester matrix" begin
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        S = sylvestermatrix(a,b)
        R = [1  4  0;
             2  5  4;
             3  0  5]
        isequal(S,R)
    end
    @test begin
        a = Polynomial([1.0,2.0,3.0],:s)
        b = Polynomial([4.0,5.0],:s)
        S = sylvestermatrix(a,b,degx=1)
        R = [1.0  0.0  4.0  0.0  0.0;
             2.0  1.0  5.0  4.0  0.0;
             3.0  2.0  0.0  5.0  4.0;
             0.0  3.0  0.0  0.0  5.0]
        isequal(S,R)
    end
end

@testset "Testing coprimeness of two polynomials" begin
    @test begin
        a = Polynomial([1,2,3],:s);
        b = Polynomial([1,2],:s);
        iscoprime(a,b)
    end
    @test begin
        a = Polynomial([1,2,3],:s);
        b = a*Polynomial([1,2],:s);
        !iscoprime(a,b)
    end
end

@testset "Testing Schur (discrete-time) stability of a polynomial" begin
    @test begin
        a = fromroots([-1/2,0, 2/3,-1/3im, 1/3im])
        isschurstable(a)
    end
    @test begin
        a = fromroots([2,0, 2/3,-3im, 3im])
        !isschurstable(a)
    end
end

@testset "Testing Hurwitz (continuous-time) stability of a polynomial" begin
    @test begin
        a = fromroots([-1, -2, -3+4im, -3-4im])
        ishurwitzstable(a)
    end
    @test begin
        a = fromroots([-1, -2, 3+4im, 3-4im])
        !ishurwitzstable(a)
    end
end

@testset "Solving the ax=b equations" begin
    @test begin
        a = Polynomial([1,2,3],:s);
        b = a*Polynomial([1,2],:s);
        x = axb(a,b)                # it should find the solution
        a*x≈b
    end
    @test begin
        a = Polynomial([1,2],:s)
        b = Polynomial([1,2,3],:s);
        x = axb(a,b)                # it should detect there is no solution
        isnothing(x)
    end
end

@testset "Solving the ax+by=0 equation" begin
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        c = Polynomial([1,1],:s)
        a = a*c
        b = b*c
        x, y = axby0(a,b)
        norm(a*x+b*y)<1e-8
    end
end

@testset "Solving the ax+by=c equation" begin
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        c = Polynomial([6,7,8],:s)
        x,y = axbyc(a,b,c)
        a*x+b*y≈c
    end
    @test begin
        a = Polynomial([1,2,3],:s)*Polynomial([1,1],:s)
        b = Polynomial([4,5],:s)*Polynomial([1,1],:s)
        c = Polynomial([5,6,7,8,9],:s)
        x,y = axbyc(a,b,c)
        a*x+b*y≈c
    end
end

@testset "Solving the ax̃+ãx=b+b̃ equation" begin
    @test begin
        a = Polynomial([-0.12, -0.29, 1],:s)
        b = Polynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:s)
        x = axxabb(a,b)
        cconj(a)*x+a*cconj(x)≈b+cconj(b)
    end
end

@testset "Solving the ax̃+b̃y=c+d̃ equation" begin
    @test begin
        a = Polynomial([-0.12, -0.29, 1],:s)
        b = Polynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:s)
        c = Polynomial([1.2, 3.4, 5.6],:s)
        d = Polynomial([3.4, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2],:s)
        x,y = axbycd(a,b,c,d)
        a*cconj(x)+cconj(b)*y≈c+cconj(d)
    end
end

@testset "Solving the equation ax+y=c while minimizing ||y||₁" begin
    @test begin
        a = Polynomial([-0.12, -0.29, 0.71, 2.79, 2.92, 1],:d)
        c = Polynomial([1.86, -0.34, -1.14, -0.21, 1.19, -1.12],:d)
        x, y, val = axycminl1(a,c)
        a*x+y≈c
    end
end
