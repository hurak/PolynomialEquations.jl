using Polynomials
using PolynomialEquations
using Test

@testset "Scaling of a polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        ρ=10
        as = scale(a,ρ)
        r = Polynomial([1,20,300,4000,50000],:s)
        isequal(as,r)
    end;
end;

@testset "Paraconjugate of a polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        p = paraconj(a)
        r = Polynomial([1,-2,3,-4,5],:s)
        isequal(p,r)
    end
    @test begin
        a = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
        p = paraconj(a)
        r = Polynomial([1-1im, -(2-2im), 3-3im, -(4-4im), 5-5im],:s)
        isequal(p,r)
    end
end;

@testset "Reciprocal polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        p = reciprocal(a)
        r = Polynomial([5,4,3,2,1],:s)
        isequal(p,r)
    end
    @test begin
        c = Polynomial([1+1im, 2+2im, 3+3im, 4+4im, 5+5im],:s)
        p = reciprocal(c)
        r = Polynomial([5-5im, 4-4im, 3-3im, 2-2im, 1-1im],:s)
        isequal(p,r)
    end
end;
