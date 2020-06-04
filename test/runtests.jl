using Laguerre
using Test

@testset "LaguerrePolynomial construction" begin
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



@testset "LaguerrePolynomial constant multiplication" begin
    orders = [0,2,3]
    coeffs = [1, 2, 2]
    L1 = LaguerrePolynomial(orders, coeffs)
    L2 = L1 * 3
    @test L2.orders == [0,2,3]
    @test L2.coeffs == 3 * coeffs
end