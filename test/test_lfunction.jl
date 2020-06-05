@testset "LaguerreFunction construction" begin
    # test construction of LaguerreFunction
    orders = [0,2,3]
    coeffs = [1.0, 2, 0.5]
    L = LaguerreFunction(orders, coeffs)
    @test L.orders == orders
    @test all(isapprox.(coeffs, L.coeffs))
    # test single Laguerre Function
    L = LaguerreFunction(2, 1 // 2)
    @test L.orders == [2]
    @test L.coeffs == [1 // 2]

    # Laguerre Function L3 with coefficient 1.
    L = LaguerreFunction(3)
    @test L.orders == [3]
    @test L.coeffs == [1]
    # test throw
    @test_throws DimensionMismatch LaguerreFunction([1,2], [2,3,2])
    @test_throws DomainError LaguerreFunction(-2, 1)
    # test string representation
    L = LaguerreFunction([2,3], [1,2])
    st = Laguerre.lashow(L)
    @test st == "(1)Lf2 + (2)Lf3"
end


@testset "LaguerreFunction addition" begin
    L1 = LaguerreFunction(2, 1 // 2)
    L2 = LaguerreFunction(3, 1.0)
    L = L1 + L2
    @test typeof(L) == LaguerreFunction{typeof(1.0)}
    @test L.orders == [2, 3]
    @test all(isapprox.(L.coeffs, [1 / 2, 1.0]))
    # 
    L1 = LaguerreFunction(2, 1)
    L2 = LaguerreFunction(2, 1 // 3)
    L = L1 + L2
    @test typeof(L) == LaguerreFunction{typeof(1 // 3)}
    @test L.orders == [2]
    @test L.coeffs == [1 + 1 // 3]
     
    # test subtraction
    L = L1 - L2
    @test L.orders == [2]
    @test L.coeffs == [1 - 1 // 3]
    
end



@testset "LaguerreFunction constant operation" begin
    orders = [0,2,3]
    coeffs = [1, 2, 2]
    L1 = LaguerreFunction(orders, coeffs)
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