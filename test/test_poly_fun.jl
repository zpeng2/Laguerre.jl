# test interoperation between LaguerrePolynomial and LaguerreFunction.
# in particular, a product of LaguerrePolynomial and LaguerreFunction 
# can be expanded into a linear combination of LaguerreFunctions.
@testset "Test throw" begin
    l1 = LaguerrePolynomial(2, 2.0)
    l2 = LaguerreFunction(3, 1.0)
    @test_throws ErrorException l1 + l2
    @test_throws ErrorException l1 - l2 
end

@testset "Product of LaguerrePolynomials" begin
    k = 5
    L2 = LaguerrePolynomial(2, 1.0)
    Lk = LaguerrePolynomial(k, 1.0)
    L = L2 * Lk
    # L2Lk = k(k-1)/2 L_{k-2} -2k(k-1)L_{k-1}
    #   +(3k^2-k)L_k -2k(k+1)L_{k+1} 
    #   (k+1)(k+2)/2 L_{k+2}
    res = k * (k - 1) / 2 * LaguerrePolynomial(k - 2)  -
    2 * k * (k - 1) * LaguerrePolynomial(k - 1) +
    (3 * k^2 - k) * LaguerrePolynomial(k) -
    2 * k * (k + 1) * LaguerrePolynomial(k + 1) +
    (k + 1) * (k + 2) / 2 * LaguerrePolynomial(k + 2)
    @test L == res
end