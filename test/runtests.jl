using PolynomialEquations
using Test

@testset "Scaling a polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        ρ=10
        as = scale(a,ρ)
        ac = Polynomial([1,20,300,4000,50000],:s)
        isequal(as,ac)
    end;
end;

end
