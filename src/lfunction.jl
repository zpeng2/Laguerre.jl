"""
A linear combination of Laguerre functions.
"""
struct LaguerreFunction{T} <: AbstractLaguerre{T}
    # representation of a linear Combination of Laguerre function
    # it's just Laguerre polynomial multiplied by exp(-x/2)
    orders::Vector{Int} # order of Laguerre polynomials.
    coeffs::Vector{T} # coefficients of Laguerre polynomials.
    len::Int # length of orders
    function LaguerreFunction(orders::Vector{Int}, coeffs::Vector{T}) where T 
        length(orders) == length(coeffs) || throw(DimensionMismatch("orders and coeffs are of different length."))
        all(orders .>= 0) || throw(DomainError("negative orders not allowed."))
        # arrange in a ascending order.
        perm = sortperm(orders)
        orders = orders[perm]
        coeffs = coeffs[perm]
        return new{T}(orders, coeffs, length(orders))
    end
end



function LaguerreFunction(order::Int, coef::Number)
    # handle a single Laguerre Function.
    return LaguerreFunction([order], [coef])
end


function LaguerreFunction(order::Int)
    # A single Laguerre Function with coefficient 1.
    return LaguerreFunction(order, 1)
end

