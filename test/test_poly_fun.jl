# test interoperation between LaguerrePolynomial and LaguerreFunction.
# in particular, a product of LaguerrePolynomial and LaguerreFunction 
# can be expanded into a linear combination of LaguerreFunctions.
@testset "Test throw" begin
    l1 = LaguerrePolynomial(2, 2.0)
    l2 = LaguerreFunction(3, 1.0)
    @test_throws ErrorException l1 + l2
    @test_throws ErrorException l1 - l2 
end

