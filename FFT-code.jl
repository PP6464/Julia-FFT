function ispow2(x::Int64)
    return x > 0 && (x & x - 1) == 0
end

function nextpow2(x::Int64)
    if ispow2(x)
        return x
    else
        return 2^(ceil(log2(x)))
    end
end

function isprime(n)
    if n <= 1
        return false
    end
    for i in 2:isqrt(n)
        if n % i == 0
            return false
        end
    end
    return true
end


# Pad a vector to a custom length with a custom padder (default is zero(T))
function pad_vector(x::Vector{T}, padder::T = zero(T), new_length::Int64, add_before::Bool = false) where T <: Any
    if add_before
        return vcat([padder for _ in 1:new_length - length(x)], x)
    else
        return vcat(x, [padder for _ in 1:new_length - length(x)])
    end
end

# Pads a matrix to custom dimensions, useful for FFT convolution 
function pad_matrix(A::Matrix{T}, m::Int64, n::Int64) where T <: Number
    A_padded = zeros(T, m, n)
    A_padded[1:size(A, 1), 1:size(A, 2)] = A
    return A_padded
end

# Returns a vector of the rows of the matrix
function rows(matrix::Matrix)
    m = size(matrix)[begin]

    return [matrix[i, begin:end] for i in 1:m]
end

# Returns a vector of the columns of the matrix
function columns(matrix::Matrix)
    n = size(matrix)[end]

    return [matrix[begin:end, i] for i in 1:n]
end

# Given a vector of rows or columns, this function will convert the vector to a matrix
# columns is used to specify whether the vector supplied is a vector of columns (true) or a vector of rows (false)
function convert_to_matrix(v::Vector{Vector{T}}, columns::Bool = false) where T <: Number
    if columns
        m, n = length(v[1]), length(v)
        matrix = zeros(T, m, n)
        for i in 1:n
            matrix[begin:end, i] = v[i]
        end
        return matrix
    else 
        m, n = length(v), length(v[1])
        matrix = zeros(T, m, n)
        for i in 1:m
            matrix[i, begin:end] = v[i]
        end
        return matrix
    end
end

function dft(x::Vector{T}) where T <: Number

end

function fft(x::Vector{T}) where T <: Number
    N = length(x)

    if ispow2(N)
        x_even = x[1:2:end]
        x_odd = x[2:2:end]

        factors = @. exp(-2pi * im * (0:N-1) / N)

        return vcat(x_even .+ x_odd .* factors[1:N÷2], x_even .+ x_odd .* factors[1+(N÷2):end])
    else if isprime(N)

    else

    end
end

function fft(x::Vector{T}) where T <: Real
    N = length(x)

end