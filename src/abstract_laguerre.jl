# abstract type for LaguerrePolynomial and LaguerreFunction
abstract type AbstractLaguerre{T <: Number} end


Base.iterate(la::AbstractLaguerre{T}, state = 1) where T = state > la.len ? nothing : ((la.orders[state], la.coeffs[state]), state + 1)


function lashow(la::AbstractLaguerre{T}) where T
    # return a string representation of a linear Combination of Laguerre polynomials
    st = ""
    if la isa LaguerrePolynomial{T} 
        name = "L"
    elseif la isa LaguerreFunction{T}
        name = "Lf"
    end
    for i in 1:la.len - 1
        st *= "($(la.coeffs[i]))" * name * "$(la.orders[i]) + "
    end
    # add the last one
    st *= "($(la.coeffs[end]))" * name * "$(la.orders[end])"
    return st
end


Base.show(io::IO, la::AbstractLaguerre{T}) where T = println(lashow(la))


# implement addition of LaguerrePolynomial or LaguerreFunction .
function Base.:+(self::AbstractLaguerre{T}, other::AbstractLaguerre{S}) where {T,S} 
    # addition of two Laguerre Polynomials or functions.
    # need to be the same type.
    nameof(typeof(self)) == nameof(typeof(other)) || throw(ErrorException("cannot add LaguerrePolynomial to a LaguerreFunction"))
    l = addla(self, other)
    return l
end


function addla(l1::AbstractLaguerre{T}, l2::AbstractLaguerre{S}) where {T,S}
    # implementation of addition of Laguerre Polynomials/Functions.
    # need to copy to create a new array.
    orders = copy(l1.orders)
    coefs = copy(l1.coeffs)
    # type promotion 
    R = promote_type(T, S)
    # convert coeffs to type R.
    coefs = R.(coefs)
    # iterate over l2 add each to l1
    for (order, coef) in l2
        if order in orders
            id = findall(orders .== order)[1]
            coefs[id] += coef
        else
            append!(orders, order)
            append!(coefs, R(coef))
        end
    end
    # get the symbol of the typename.
    typename = nameof(typeof(l1))
    # constructor expression.
    expr = :($typename($orders, $coefs))
    return eval(expr)
end



# multiplication of LaguerrePolynomial/LaguerreFunction with a number.
function Base.:*(l::AbstractLaguerre{T}, num::S) where {S <: Number} where T
    # product of LaguerrePolynomial with a constant.
    R = promote_type(T, S)
    coefs = R.(copy(l.coeffs))
    coefs *= num
    # get the symbol of the typename.
    typename = nameof(typeof(l))
    # constructor expression.
    expr = :($typename($(l.orders), $coefs))
    return eval(expr)
end


Base.:*(num::S, l::AbstractLaguerre{T}) where T where S <: Number = l * num


function Base.:-(l1::AbstractLaguerre{T}, l2::AbstractLaguerre{S}) where {T,S}
    # subtraction of two LaguerrePolynomial or LaguerreFunction
    # need to be the same type.
    # no need to check it here, because addition will check .
    R = promote_type(T, S)
    minus1 = R(-1)
    return l1 + minus1 * l2
end




function Base.:/(la::AbstractLaguerre{T}, num::S) where {S <: Number} where T
    num != 0 || throw(ErrorException("division by zero not allowed."))
    # divide a LaguerrePolynomial by a number can be done.
    R = promote_type(T, S)
    coefs = R.(copy(la.coeffs))
    coefs /= num
    # get the symbol of the typename.
    typename = nameof(typeof(la))
    # constructor expression.
    expr = :($typename($(la.orders), $coefs))
    return eval(expr)
end


# check equality of two LaguerrePolynomials or LaguerreFunctions
function Base.:(==)(l1::AbstractLaguerre{T}, l2::AbstractLaguerre{S}) where {T,S}
    if (nameof(typeof(l1)) != nameof(typeof(l2)))
        return false
    end
    return (l1.orders == l2.orders) && (l1.coeffs == l2.coeffs)
end



function simplify(l::AbstractLaguerre{T}) where T
    # remove zero elements from LaguerrePolynomial or LaguerreFunction
    # find nonzero coeffs
    idx = findall(l.coeffs .!= 0)
    if isempty(idx)
        # all zeros.
        # just return the zero order term with 0 coefficient.
        coeffs = 0
        orders = 0
    elseif length(idx) == l.len
        # all coefficients are nonzero.
        # nothing needs to be done.
        return l
    else
        # some are zero.
        coeffs = l.coeffs[idx]
        orders = l.orders[idx]
    end
    # get the symbol of the typename.
    typename = nameof(typeof(l))
    # constructor expression.
    expr = :($typename($orders, $coeffs))
    return eval(expr)
end