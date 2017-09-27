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


"""
Interpolation module for polynomial interpolation real valued data.
"""
module Interpolation

using Main: unit_pulse, shift

export polynomial_interpolate


"""
    polynomial_interpolate(x::Array{<:Real, 1}, y::Array{<:Real, 1})

Polynomial interpolation of `x` and `y` vectors. Degree of the interpolated
polynomial is `n - 1` where `n` is the length of `x`.
"""
function _polynomial_interpolate(x::Array{<:Real, 1},
                                 y::Array{<:Real, 1})
    n = length(x)
    deg = n - 1

    mat = zeros(n, n)
    for i = 0 : deg
        mat[:, i + 1] = x.^i
    end
    a = zeros(n)
    try
        a = inv(mat) * y
    catch
        throw(ErrorException("No unique solution for the intepolaton"))
    end

    return t -> dot(a, [t^i for i = 0 : deg])

end  # End of _polynomial_interpolate function


"""
    polynomial_interpolate(x::Array{<:Real, 1},
                           y::Array{<:Real, 1};
                           deg::Union{Int, Bool}=nothing)

Polynomial interpolation of real valued points `x` and `y`. The degree of  the
interpolation polynomial may be any integer less than the length `n` of `x`.
"""
function polynomial_interpolate(x::Array{<:Real, 1},
                                y::Array{<:Real, 1};
                                deg::Union{Int, Bool}=nothing)

    if ~(length(x) == length(y))
        throw(ArgumentError("Vector lengths does not match"))
    end

    # If the interpolation points are mixed, sort them.
    xy = sortrows([x y])
    x, y = xy[:, 1], xy[:, 2]

    n = length(x)
    if deg == nothing || deg == n - 1
        f_interpolated = _polynomial_interpolate(x, y)
        function func_full(t)
            if t < x[1] || t > x[end]
                throw(DomainError("$t is out of interpolaton domain $(x[1]), $(x[end])"))
            elseif t == x[1]
                return y[1]
            elseif t == x[end]
                return y[end]
            else
                return f_interpolated(t)
            end  # End of if-elseif-else
        end  # End of func function
        return func_full
    elseif 1 <= deg < n - 1
        if mod(n - 1, deg) == 0
            num_intevals = floor(Int, (n - 1) / deg)
        else
            num_intevals = floor(Int, (n - 1) / deg) + 1
        end
        funcs = Array{Function}(num_intevals)
        mask = falses(n)
        mask[1 : deg + 1] = true
        for k = 1 : num_intevals
            x_block = x[mask]
            y_block = y[mask]
            # Pad the intervals whose lengths are not equal to deg + 1
            # This padding is done by adding a perturbation (1e-6 here) to
            #   the last element in the interval.
            if ~(length(x_block) == deg + 1)
                x_pad = [x_block[end] + 1e-6 for i = 1 : deg + 1 - length(x_block)]
                y_pad = [y_block[end] + 1e-6 for i = 1 : deg + 1 - length(y_block)]
                x_block = vcat(x_block, x_pad)
                y_block = vcat(y_block, y_pad)
            end  # End of if
            fk = _polynomial_interpolate(x_block, y_block)
            funcs[k] = t -> fk(t) * unit_pulse(t, tl=x_block[1], th=x_block[end])
            mask = shift(mask, deg)
        end  # End of for

        function func_piecewise(t)
            if t < x[1] || t > x[end]
                throw(DomainError("$t is out of interpolation domain $(x[1]), $(x[end])"))
            elseif t == x[1]
                return y[1]
            elseif t == x[end]
                return y[end]
            else
                return sum([funcs[k](t) for k = 1 : num_intevals])
            end  # End of if-elseif-else
        end  # End of func
        return func_piecewise
    else
        throw(ArgumentError("For $n points expected degree is less than $n, got $deg instead"))
    end  # End of if-elseif-else.

end  # End of polynomial_interpolate function

end  # End of the Interpolation module
