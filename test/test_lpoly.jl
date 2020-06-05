@testset "LaguerrePolynomial" begin
    # test construction of LaguerrePolynomial
    orders = [0,2,3]
    coeffs = [1.0, 2, 0.5]
    L = LaguerrePolynomial(orders, coeffs)
    @test L.orders == orders
    @test all(isapprox.(coeffs, L.coeffs))
    # test single Laguerre Polynomial
    L = LaguerrePolynomial(2, 1 // 2)
    @test L.orders == [2]
    @test L.coeffs == [1 // 2]

    # Laguerre Polynomial L3 with coefficient 1.
    L = LaguerrePolynomial(3)
    @test L.orders == [3]
    @test L.coeffs == [1]
    # test throw
    @test_throws DimensionMismatch LaguerrePolynomial([1,2], [2,3,2])
    @test_throws DomainError LaguerrePolynomial(-2, 1)
    # test string representation
    L = LaguerrePolynomial([2,3], [1,2])
    st = Laguerre.lashow(L)
    @test st == "(1)L2 + (2)L3"
end


@testset "LaguerrePolynomial addition" begin
    L1 = LaguerrePolynomial(2, 1 // 2)
    L2 = LaguerrePolynomial(3, 1.0)
    L = L1 + L2
    @test typeof(L) == LaguerrePolynomial{typeof(1.0)}
    @test L.orders == [2, 3]
    @test all(isapprox.(L.coeffs, [1 / 2, 1.0]))
    # 
    L1 = LaguerrePolynomial(2, 1)
    L2 = LaguerrePolynomial(2, 1 // 3)
    L = L1 + L2
    @test typeof(L) == LaguerrePolynomial{typeof(1 // 3)}
    @test L.orders == [2]
    @test L.coeffs == [1 + 1 // 3]
     
    # test subtraction
    L = L1 - L2
    @test L.orders == [2]
    @test L.coeffs == [1 - 1 // 3]
    
end



@testset "LaguerrePolynomial * number" begin
    orders = [0,2,3]
    coeffs = [1, 2, 2]
    L1 = LaguerrePolynomial(orders, coeffs)
    L2 = L1 * 3
    # multiplication by a number
    @test L2.orders == [0,2,3]
    @test L2.coeffs == 3 * coeffs
    # division by a number
    @test_throws ErrorException L1 / 0
    L3 = L1 / 2
    @test L3.orders == L1.orders
    @test L3.coeffs == L1.coeffs / 2
end

@testset "LaguerrePolynomials product" begin
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

    # L1Lk = k L_{k-1} -2k L_k +(k+1) L_{k+1}
    L1 = LaguerrePolynomial(1, 1.0)
    Lk = LaguerrePolynomial(k, 1.0)
    L = L1 * Lk
    res = k * LaguerrePolynomial(k - 1) - 2 * k * LaguerrePolynomial(k) + (k + 1) * LaguerrePolynomial(k + 1)
    @test L ==  res
    @test Lk * LaguerrePolynomial(0) == Lk
end