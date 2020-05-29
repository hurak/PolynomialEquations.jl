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
        r = LaurentPolynomial(5:-1:0,-4:0,:z)
        isequal(p,r)
    end
    @test begin
        a = = LaurentPolynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:z)
        d = dconj(a)
        r = LaurentPolynomial([5-5im, 4-4im, 3-3im, 2-2im, 1-1im],-4:0,:z)
        isequal(p,r)
    end
end

@testset "Building a lower triangular banded Toeplitz matrix" begin
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
end

@testset "Building a Sylvester matrix" begin
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        S = sylvestermatrix(a,b)
        R = [1.0  4.0  0.0;
             2.0  5.0  4.0;
             3.0  0.0  5.0]
        isequal(S,R)
    end
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        S = sylvestermatrix(a,b,degx=1)
        R = [1.0  0.0  4.0  0.0  0.0;
             2.0  1.0  5.0  4.0  0.0;
             3.0  2.0  0.0  5.0  4.0;
             0.0  3.0  0.0  0.0  5.0]
        isequal(S,R)
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
        norm(a*x+b*y) < 1e-8
    end
end

@testset "Solving the ax+by=c equation" begin
    @test begin
        a = Polynomial([1,2,3],:s)
        b = Polynomial([4,5],:s)
        c = Polynomial([6,7,8],:s)
        x,y = axbyc(a,b,c)
        norm(a*x+b*y-c) < 1e-8
    end
end
