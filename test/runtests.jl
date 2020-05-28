using PolynomialEquations
using Test

@testset "Scaling a polynomial" begin
    @test begin
        a = Polynomial(1:5,:s)
        ρ=10
        as = scale(a,ρ)
        Polynomial(1 + 20*s + 300*s^2 + 4000*s^3 + 50000*s^4)
        isequal(a,as)
    end;
end;

end
