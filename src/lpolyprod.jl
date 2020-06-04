
function laprodexpand(r::Int, s::Int)
    # expand the product of two Laguerre polynomials into a linear
    # combination of Laguerre polynomials.
    # Lr(x)Ls(x) into sum_{t}C_{rst}Lt(x)
    orders = abs(r - s):(r + s) |> collect
    coeffs = [Crst(r, s, t) for t in orders]
    return orders, coeffs
end


function Crst_direct(r::Int, s::Int, t::Int)
    # direct method to calculate Crst
    # see Products of Laguerre Polynomials
    # Joseph Gillis and  George Weiss
    # Mathematics of Computation, Vol. 14, No. 69 (Jan., 1960), pp. 60-63
    # this function implements equation (7)
    # this is computationally expensive.
    # use recurrence relation instead.
    (r >= 0 && s >= 0 && t >= 0) || throw("negative orders not allowed.")
    # t needs to be in |r-s| and r+s
    # otherwise the coef is zero
    if t < abs(r - s) || t > r + s
        return 0
    end

    p = r + s - t
    # define nmax and nmin and find them
    nmax = r
    if s < nmax
        nmax =  s
    end
    if p < nmax
        nmax = p
    end
    nmin = ceil(Int, p / 2)
    # use rational numbers for now to keep it exact.
    C = 0 // 1
    for n in nmin:nmax
        C += (2 // 1)^(2 * n) * factorial(r + s - n) / ( factorial(r - n) * factorial(s - n) * factorial(2n - p) * factorial(p - n))
    end
    # prefactor
    C *= (-1 // 2)^p
end


function Crst(r::Int, s::Int, t::Int)
    # this implements equation (22) and (23) in the same paper.
    if r < 0 || s < 0 || t < 0 
        return 0 // 1
    end
    # (r >= 0 && s >= 0 && t >= 0) || throw("negative orders not allowed.")
    # notice that C_{rst} is symmetric with respect to any swapping of indices.
    # take r to be the larer one of r and s
    # the recurrence is given for t+1
    # shift it to calculate for t;
    if r < s
        r, s = s, r # swap.
    end
    # t needs to be in |r-s| and r+s
    # otherwise the coef is zero
    if t < r - s || t > r + s
        return 0 // 1
    end

    sigma = r + s + 1 // 1
    delta = r - s
    if t == delta
        return Rational(binomial(r, s))
    end

    # handle general case.
    f1 = t * (delta^2 - t^2)
    f2 = (2t - 1) * (delta^2 - (t - 1)^2 + t * (sigma - t ))
    f3 = -(t - 1) * (delta^2 - (t - 2)^2 + (4t - 3) * (sigma - t + 1))
    f4 = 2(t - 1) * (t - 2) * (sigma - t + 2)
    return (f2 * Crst(r, s, t - 1) + f3 * Crst(r, s, t - 2) + f4 * Crst(r, s, t - 3) ) / f1
end




# function Base.:*(l1::La, l2::La)
#     # take two linear combinations of Laguerre polynomials.
#     # and reexpand into a single Laguerre polynomial series.
#     lst = []

#     for (r, c1) in l1
#         for (s, c2) in l2
#             # construct the product 
#             orders, coeffs = laprodexpand(r, s)
#             # convert to float
#             coeffs = Float64.(coeffs)
#             # add this to the list
#             la = La(orders, coeffs)
#             la = la * c1 * c2
#             push!(lst, la)
#         end
#     end
#     # add elements in the list.
#     la = lst[1]
#     if length(lst) > 1
#         for i = 2:length(lst)
#             la += lst[i]
#         end
#     end
#     return la
# end