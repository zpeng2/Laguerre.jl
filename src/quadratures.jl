function LGnodes(alpha::Real, N::Int)
    # generate nodes of Laguerre-Gauss quadrature nodes.
    alpha > -1 || error("alpha <= -1 not allowed.")
    j = 1:N |> collect
    b = j .* (j .+ alpha)
    # +- diagonal
    dl = -sqrt.(b) 
    j = 0:N |> collect
    # main diagonal
    d = 2.0 * j .+ alpha .+ 1
    A = SymTridiagonal(d, dl)
    return eigvals(A)
end



function LGRnodes(N::Int)
    # Laguerre-Gauss-Radau quadrature nodes.
    xi = LGnodes(1.0, N - 1)
    x = [0; xi]
    return x
end


function LGR(N::Int)
    # get nodes and weights for LGR.
    x = LGRnodes(N)
    LN = eval_laguerre_function.(N, x)
    # equation 7.37 in Shen, Tang & Wang.
    w = 1 ./ ( (N + 1) * LN.^2)
    return x, w
end



function eval_laguerre_function(n::Int, x::Real)
    # evaluate Laguerre function using the recurrence relation
    L0 = exp(-x / 2)
    L1 = (1 - x) * exp(-x / 2)
    if n == 0
        return L0
    end
    if n == 1
        return L1
    end
    # n >= 2 case
    Lnew = 0.0
    for i = 2:n
        Lnew = (2 - 1 / i - x / i) * L1 - (1 - 1 / i) * L0
        L0 = L1
        L1 = Lnew 
    end
    return Lnew
end

function eval_LaguerreFunction(n::Int, x::Real)
    # method in Shen, SIAM J NUMER ANAL, Vol 38, No 4, pp. 1113-1133
    # STABLE AND EFFICIENT SPECTRAL METHODS IN UNBOUNDED
    # DOMAINS USING LAGUERRE FUNCTIONS
    signL0 = 1
    d0 = 1
    signL1 = sign(1 - x)
    d1 = abs(1 - x)
    if n == 0 
        return exp(-x / 2)
    end
    if n == 1
        return signL1 * d1
    end
    # TODO


end

function LGRquad(f::Function, N::Int)
    # numerical integration using Laguerre-Gauss-Radau quadrature.
    # from 0 to infinity.
    x, w = LGR(N)
    return sum(f.(x) .* w)
end


function laguerre_transform(f::Function, N::Int)
    # interpolation: represent f in a series expansion
    # f = sum_{n=0}^N a_n Lnhat(x)
    # where Lnhat is the Laguerre function.
    # an = \sum_{j=0}^N f(x_j)Lnh(x_j)what_j
    x, w = LGR(N)
    an = zeros(N + 1)
    for i = 1:N + 1
        Ln = eval_laguerre_function.(i - 1, x)
        an[i] = sum(f.(x) .* Ln .* w)
    end
    return an
end


function inverse_laguerre_transform(an::AbstractVector, x::AbstractArray)
    # take an array of coefficients and reconstruct the function.
    N = length(an) - 1
    # f = sum_{n=0}^N an Lnh(x)
    f = zeros(size(x))

    for (a, n) in zip(an, 0:N)
        Lnh = eval_laguerre_function.(n, x)
        f += a * Lnh
    end
    return f
end


