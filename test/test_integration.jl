@testset "Laguerre integration" begin
    # forward transform
    # this is just Lnh for n=2.
    # need to get a2 = 1 and other coeffs 0.
    f(x) = (x^2 - 4 * x + 2) * exp(-x / 2) / 2
    N = 3
    an = laguerre_transform(f, N)
    anexact = [0; 0; 1; 0]
    @test all(isapprox.(an, anexact, atol = 1e-8))

    # now inverse transform
    x = 0:0.1:10 |> collect
    fv = inverse_laguerre_transform(an, x)
    @test norm(fv - f.(x)) < 1e-10

    # test integration using Laguerre-Gauss-Radau quadrature.
    f1(x) = sin(x)^2 * exp(-x)  # -> 2/5
    N = 30
    @test isapprox(LGRquad(f1, N), 2 / 5, atol = 1e-3, rtol = 1e-4)
end