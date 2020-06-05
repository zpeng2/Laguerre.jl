
"""
A linear combination of Laguerre polynomials.
"""
struct LaguerrePolynomial{T} <: AbstractLaguerre{T}
    # representation of a linear Combination of Laguerre polynomials
    orders::Vector{Int} # order of Laguerre polynomials.
    coeffs::Vector{T} # coefficients of Laguerre polynomials.
    len::Int # length of orders
    function LaguerrePolynomial(orders::Vector{Int}, coeffs::Vector{T}) where T 
        length(orders) == length(coeffs) || throw(DimensionMismatch("orders and coeffs are of different length."))
        all(orders .>= 0) || throw(DomainError("negative orders not allowed."))
        # arrange in a ascending order.
        perm = sortperm(orders)
        orders = orders[perm]
        coeffs = coeffs[perm]
        return new{T}(orders, coeffs, length(orders))
    end
end


function LaguerrePolynomial(order::Int, coef::Number) 
    # handle a single Laguerre polynomial.
    return LaguerrePolynomial([order], [coef])
end


function LaguerrePolynomial(order::Int)
    # A single Laguerre Polynomial with coefficient 1.
    return LaguerrePolynomial(order, 1)
end



# can add LaguerrePolynomial to a constant, because a constant is just const*L0
function Base.:+(la::LaguerrePolynomial{T}, num::S) where {S <: Number} where T
    R = promote_type(T, S)
    L0 = LaguerrePolynomial(0, num)
    return la + L0
end

# order doesn't matter in addition.
function Base.:+(num::S, la::LaguerrePolynomial{T}) where T where {S <: Number} 
    return la + num
end





