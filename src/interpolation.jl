"""
Interpolation module for polynomial interpolation real valued data.
"""
module Interpolation

include("./utilities.jl")


export polynomial_interpolate, spline_interpolate, fractal_interpolate


"""
    polynomial_interpolate(x::Array{Real, 1}, y::Array{<:Real, 1})

Polynomial interpolation of `x` and `y` vectors. Degree of the interpolated
polynomial is `n - 1` where `n` is the length of `x`.
"""
function _polynomial_interpolate(x::Vector{<:Real},
                                 y::Vector{<:Real})
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
    polynomial_interpolate(x::Array{Real, 1},
                           y::Array{Real, 1};
                           deg::Union{Int, Bool}=nothing)

Polynomial interpolation of real valued points `x` and `y`. The degree of  the
interpolation polynomial may be any integer less than the length `n` of `x`.
"""
function polynomial_interpolate(x::Vector{<:Real},
                                y::Vector{<:Real};
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


"""
    spline_interpolate(x::Vector{Real}, y::Vector{Real}; spline_type::String="linear")

Spline interpolation of real vectors `x` and `y`.

Interpolation is performed by fitting of polynomial between each pairs of
`x[i] : x[i +1]`. `spline_type` determines the degree of the polynomial.
Possible values for `spline_type` is

    * `linear` : A linear polynomial `p(x) = a + bx` is interpolated.
    * `quadratic` : A quadratic polynomial `p(x) = a + bx + x^2` is interpolated.
    * `qubic` : A qubic polynomial `p(x) = a + bx + cx^2 + dx^3` is interpolated.
"""
function spline_interpolate(x::Vector{<:Real}, y::Vector{<:Real}; spline_type::AbstractString="linear")
    # Check the data length
    if ~(length(x) == length(y))
        throw(ArgumentError("Vector lengths does not match"))
    end

    # If the interpolation points are mixed, sort them.
    xy = sortrows([x y])
    x, y = xy[:, 1], xy[:, 2]

    # Define the mask function.
    # The mask function is used for boolean indexing of splines.
    function mask(t)
        ind = Vector{Bool}(n)
        for k = 1 : n
            ind[k] = Bool(unit_pulse(t, tl=x[k], th=x[k + 1]))
        end
        return ind
    end

    if spline_type == "linear"
        # Compute spline coefficients
        n = length(x) - 1
        delta_x = diff(x)
        delta_y = diff(y)
        b = \(diagm(delta_x), delta_y)
        a = y[2:end] - b .* x[2:end]

        # Construct splines
        splines = Vector{Function}(n)
        for k = 1 : n
            splines[k] = t -> a[k] + b[k] * t
        end

    elseif spline_type == "quadratic"
        # 1st derivative vector z at interpolation points x
        n = length(x) - 1
        mat = full(Tridiagonal(ones(n), ones(n + 1), zeros(n)))
        d = vcat([0], 2 * diff(y) ./ diff(x))
        println("mat: $mat d: $d")
        z = \(mat, d)

        # Compute the spline coefficients
        delta_x = diff(x)
        delta_y = diff(y)
        delta_z = diff(z)
        x_bar = x[1 : end - 1]
        y_bar = y[1 : end - 1]
        z_bar = z[1 : end - 1]
        a = -z_bar .* x_bar + (delta_z ./ (2 * delta_x)) .* (x_bar.^2) + y_bar
        b = z_bar - (delta_z ./ delta_x) .* x_bar
        c = delta_z ./ (2 * delta_x)

        # Construct splines
        splines = Vector{Function}(n)
        for k = 1 : n
            splines[k] = t -> a[k] + b[k] * t + c[k] * t^2
        end

    elseif spline_type == "cubic"
        # 2nd derivative vector z at interpolation points x
        n = length(x) - 1
        delta_x = diff(x)
        delta_y = diff(y)
        di = 2 * ([delta_x[1]; delta_x] + [0; delta_x[2 : end]; 0])
        mat = full(Tridiagonal(delta_x, di, delta_x))
        d = 6 * [delta_y[1] / delta_x[1];
                delta_y[2 : end] ./ delta_x[2 : end] - delta_y[1 : end - 1] ./ delta_x[1 : end - 1];
                -delta_y[end] / delta_x[end]]
        z = \(mat, d)

        # Compute spline coefficients
        x_bar = x[1 : end - 1]
        y_bar = y[1 : end - 1]
        z_bar = z[1 : end - 1]
        x_tilde = x[2 : end]
        y_tilde = y[2 : end]
        z_tilde = z[2 : end]
        k1 = z_bar ./ (6 * delta_x)
        k2 = z_tilde ./ (6 * delta_x)
        k3 = (y_tilde ./ delta_x) - (z_tilde .* delta_x / 6)
        k4 = (y_bar ./ delta_x) - (z_bar .* delta_x / 6)
        a = k1 .* (x_tilde.^3) - k2 .* (x_bar.^3) - k3 .* x_bar + k4 .* x_tilde
        b = -3 * k1 .* (x_tilde.^2) + 3 * k2 .* (x_bar.^2) + k3 - k4
        c = 3 * k1 .* x_tilde - 3 * k2 .* x_bar
        d = -k1 + k2

        # Construct splines
        splines = Vector{Function}(n)
        for k = 0 : n - 1
            splines[k + 1] = t -> a[k + 1] + b[k + 1] * t + c[k + 1] * t^2 + d[k + 1] * t^3
        end
    end  # End of if-elseif-else.

    # Constuct the interpolation polynomial i.e. splines.
    function func(t)
        if t < x[1] || t > x[end]
            throw(DomainError("$t is out of interpolation domain $(x[1]), $(x[end])"))
        elseif t == x[end]
            return y[end]
        else
            return splines[mask(t)][1](t)
        end  # End of if-elseif-else
    end  # End of func function.

    return func
end  # End of spline_interpolate function.


"""
function fractal_interpolate(x::Vector{Real},
                             y::Vector{Real},
                             d::Vector{Real},
                             func0::Function;
                             num_iter::Integer=10)

Fractal interpolation of real points `x` and `y`. `d` determines the traslation
and has to have elements between 0 and 1. `func0` is the initial function to
interpolate. `num_iter` is the number of iterations to have calculate the
fractal interpolation function.
"""
function fractal_interpolate(x::Vector{<:Real},
                             y::Vector{<:Real},
                             d::Vector{<:Real},
                             func0::Function;
                             tol::AbstractFloat=1e-3,
                             num_iter::Int=10)
    # Check the data length.
    if length(x) != length(y)
        throw(ArgumentError("Vector lengths does not match"))
    end

    # Check translation vector d.
    if !all(0 .< d .< 1)
        error("Elements of d must be between 0 and 1.")
    end

    # If the interpolation points are mixed, sort them.
    xy = sortrows([x y])
    x, y = xy[:, 1], xy[:, 2]

    # Compute the transformation coefficients
    n = length(x) - 1
    delta_x = diff(x)
    delta_y = diff(y)
    a = delta_x / (x[end] - x[1])
    c = delta_y / (x[end] - x[1]) - d * (y[end] - y[1]) / (x[end] - x[1])
    e = x[end] / (x[end] - x[1]) * x[1 : end - 1] - x[1] / (x[end] - x[1]) * x[2 : end]
    f = x[end] / (x[end] - x[1]) * y[1 : end - 1] - x[1] / (x[end] - x[1]) * y[2 : end] -
        d * (x[end] * y[1] - x[1] * y[end]) / (x[end] - x[1])

    # Define functional transformation.
    function tf(func)
        indexes = 1 : length(x)
        function ff(t)
            index = t .> x
            if all(index)
                return y[end]
            elseif any(index)
                k = indexes[index][end]
                return c[k] * ((t - e[k]) / a[k]) + d[k] * func((t - e[k]) / a[k]) + f[k]
            else
                return y[1]
            end  # if-elseif-else
        end  # ff
        return ff
    end  # tf

    # Calculate number of iterations
    sigma = maximum(d)
    func1 = tf(func0)
    diff_func(x) = func1(x) - func0(x)
    # dist = maximum(abs(map(func1, x) - map(func0, x)))
    dist = fnormp(diff_func, x[1], x[end], 2)
    num_iter = ceil(Int, log(tol * (1 - sigma) / dist) / log(sigma))

    # Iterate the initial function.
    func = func0
    for i = 1 : num_iter
        func = tf(func)
    end
    return func
end  # fractal_interpolate


function hidden_variable_fractal_interpolate(x::Vector{<:Real},
                                             y::Vector{<:Real},
                                             d::Vector{<:Real},
                                             h::Vector{<:Real},
                                             l::Vector{<:Real},
                                             m::Vector{<:Real},
                                             hh::Vector{<:Real},
                                             func0::Function;
                                             tol::AbstractFloat=1e-3,
                                             num_iter::Int=10)

    # Check the data length.
    if length(x) != length(y)
        throw(ArgumentError("Vector lengths does not match"))
    end

    # Check translation vector d.
    # TODO: Check the hidden variables are less than.
    if !all((0 .< d .< 1) .& (0 .< h .< 1) .& (0 .< l .< 1) .& (0 .< m .< 1))
        error("Elements of hidden variables must be between 0 and 1.")
    end

    # If the interpolation points are mixed, sort them.
    xy = sortrows([x y])
    x, y = xy[:, 1], xy[:, 2]

    # Compute the transformation coefficients
    p = y[1 : end - 1] - d * y[1] - h * hh[1]
    q = h[1 : end - 1] - l * y[1] - m * hh[1]
    r = y[2 : end] - d * y[end] - h * hh[end]
    s = h[2 : end] - l * y[end] - m * hh[end]

    a = (x[2 : end] - x[1 : end - 1]) / b
    e = ((x[end] .* x[1 : end - 1])  - (x[1] - x[2 : end])) / b
    c = (r - p) / b
    k = (s - q) / b
    f = p - c * x[1]
    g = q - k * x[1]


    # Define functional transformation.
    function tf(func)
        indexes = 1 : length(x)
        function ff(t)
            index = t .> x
            if all(index)
                return y[end]
            elseif any(index)
                k = indexes[index][end]
                return c[k] * ((t - e[k]) / a[k]) + d[k] * func((t - e[k]) / a[k]) + f[k]
            else
                return y[1]
            end  # if-elseif-else
        end  # ff
        return ff
    end  # tf

    #  # Calculate number of iterations
    #  sigma = maximum(d)
    #  func1 = tf(func0)
    #  diff_func(x) = func1(x) - func0(x)
    #  # dist = maximum(abs(map(func1, x) - map(func0, x)))
    #  dist = fnormp(diff_func, x[1], x[end], 2)
    #  num_iter = ceil(Int, log(tol * (1 - sigma) / dist) / log(sigma))

    # Iterate the initial function.
    func = func0
    for i = 1 : num_iter
        func = tf(func)
    end
    return func

end  # hidden_variable_fractal_interpolate

end  # End of the Interpolation module
