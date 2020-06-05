
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


function LaguerrePolynomial(order::Int, coef::T) where T 
    # handle a single Laguerre polynomial.
    return LaguerrePolynomial([order], [coef])
end


function LaguerrePolynomial(order::Int)
    # A single Laguerre Polynomial with coefficient 1.
    return LaguerrePolynomial(order, 1)
end




# function Base.:+(self::LaguerrePolynomial{T}, other::LaguerrePolynomial{S}) where {T,S} 
#     # addition of two Laguerre Polynomials.
#     l = addla(self, other)
#     return l
# end


# function addla(l1::LaguerrePolynomial{T}, l2::LaguerrePolynomial{S}) where {T,S}
#     # implementation of addition of Laguerre Polynomials
#     # need to copy to create a new array.
#     orders = copy(l1.orders)
#     coefs = copy(l1.coeffs)
#     # type promotion 
#     R = promote_type(T, S)
#     # convert coeffs to type R.
#     coefs = R.(coefs)
#     # iterate over l2 add each to l1
#     for (order, coef) in l2
#         if order in orders
#             id = findall(orders .== order)[1]
#             coefs[id] += coef
#         else
#             append!(orders, order)
#             append!(coefs, R(coef))
#         end
#     end
#     return LaguerrePolynomial(orders, coefs)
# end


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


# function Base.:*(l::LaguerrePolynomial{T}, num::S) where {S <: Number} where T
#     # product of LaguerrePolynomial with a constant.
#     R = promote_type(T, S)
#     coefs = R.(copy(l.coeffs))
#     coefs *= num
#     return LaguerrePolynomial(l.orders, coefs)
# end

# order doesnt matter.
# Base.:*(num::S, l::LaguerrePolynomial{T}) where T where S <: Number = l * num


# function Base.:-(l1::LaguerrePolynomial{T}, l2::LaguerrePolynomial{S}) where {T,S}
#     # subtraction of two Laguerre Polynomials.
#     R = promote_type(T, S)
#     minus1 = R(-1)
#     return l1 + minus1 * l2
# end



# function Base.:/(la::LaguerrePolynomial{T}, num::S) where {S <: Number} where T
#     num != 0 || throw(ErrorException("division by zero not allowed."))
#     # divide a LaguerrePolynomial by a number can be done.
#     R = promote_type(T, S)
#     coefs = R.(copy(la.coeffs))
#     coefs /= num
#     return LaguerrePolynomial(la.orders, coefs)
# end


