# This file defines some utility fuctions

using QuadGK

"""
    unit_step(t::Real)

Unit step function ``u(t) = 0.5 * sgn(t) + 0.5```
"""
function unit_step(t::Real)
    if t == 0
        return 1
    end
    return 1 / 2 * sign(t) + 1 / 2
end  # End of unit_step function


"""
    unit_pulse(t::Real; tl::Real=0, th::Real=1)

Unit pulse function with pulse length (th - t1)
"""
function unit_pulse(t::Real; tl::Real=0, th::Real=1)
    if tl >= th
        throw(ArgumentError("tl must be less than th."))
    end
    return unit_step(t - tl) - unit_step(t - th)
end  # End of unit_pulse function


"""
    shift(vec::Array{<:Real, 1}, k::Union{Bool, Number})

Shift array `vec` for `k` positions to the right and
append zeros or falses from the left.
"""
function shift(vec::BitArray{1}, k::Int)
    temp = vec[1 : end - k]
    if typeof(vec[1]) == Bool
        return [falses(k); temp]
    else
        return [zeros(k); temp]
    end
end  # End of shift function


function fnorm1(f::Function, t0::Real, tf::Real)
    integrand(x) = abs(f(x))
    return quadgk(integrand, t0, tf)[1]
end


function fnormInf(f::Function, t0::Real, tf::Real)
    domain = linspace(t0, tf, 1000)
    return maximum(f.(domain))
end


function fnormp(f::Function, t0::Real, tf::Real, p::Real=1)
    if p == 1
        return fnorm1(f, t0, tf)
    elseif p == Inf
        return fnormInf(f, t0, tf)
    else
        integrand(x) = abs(f(x))^p
        return quadgk(integrand, t0, tf)[1]^(1 / p)
    end
end
