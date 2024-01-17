using Primes

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

# Pad a vector to a custom length with a custom padder (default is zero(T))
function pad_vector(x::Vector{T}, new_length::Int64, padder::T = zero(T), add_before::Bool = false) where T <: Any
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

# Get prime factors of a positive integer
function primefactors(n::Integer)
    if n <= 0
        error("Argument 'n' to be factored must be a positive integer.")
    end

    f = factor(n)
    return sort(collect(keys(f)))
end

# Find the smallest primitive root modulo m
function primroot(m::Integer)
    if !isprime(m)
        error("Argument 'm' must be a prime number")
    end
    if m == 2; return 1; end

    P = primefactors(m-1)
    for r = 2:(m-1)
        not_found = true
        for p in P
            if powermod(r, div(m-1, p), m) == 1
                not_found = false
            end
        end
        if not_found
            return r
        end
    end

    return 0
end

function primitive_roots_of(N::Int64)
    g1 = primroot(N)
    g2 = other_primitive_root(g1, N)

    return g1, g2
end

# Get the other primitive root modulo m (multiplicative inverse of the smaller root)
function other_primitive_root(smaller_root::Int64, m::Int64)
    # Only being used in Rader FFT so m is prime
    if !isprime(m)
        throw(ArgumentError("m must be a prime number"))
    end

    return (smaller_root^(m - 2)) % m
end

function rader_fft(x::Vector{T}) where T <: Number
    if !isprime(x)
        throw(ArgumentError("Length for Rader FFT input should be a prime number"))
    end

    N = length(x)
    g1, g2 = primitive_roots_of(N)

    ω = exp(-2pi * im / N)

    result = zeros(Complex{Float64}, N)
    # We know the first value in the FFT output must simply be the sum of the entries of x
    result[1] = sum(x)

    # Construct the DFT matrix
    dft_matrix = zeros(Complex{Float64}, N-1, N-1)

    for i in 1:N
        for j in 1:N
            dft_matrix[i,j] = ω^(ij % N)
        end
    end

    dft_matrix = dft_matrix[2:end, 2:end]

    # Construct the permutation matrices
    # Get the pos of 1 in each column
    g1_pos_ones = [g1^i for i in 0:N-1]
    g2_pos_ones = [g2^i for i in 0:N-1]

    # Construct the columns using the information
    g1_perm_columns = [[i == pos_one ? 1 : 0  for i in 0:N-1] for pos_one in g1_pos_ones]
    g2_perm_columns = [[i == pos_one ? 1 : 0 for i in 0:N-1] for pos_one in g2_pos_ones]

    # Construct the permutation matrices
    g1_perm = convert_to_matrix(g1_perm_columns, true)
    g2_perm = convert_to_matrix(g2_perm_columns, true)

    # Make a circular permutated version of DFT matrix to make computations faster
    dft_circular = g1_perm * dft_matrix * g2_perm

    # Find order of outputs
    y_indices = [i for i in 1:N]
    y_indices_perm1 = []    
end

function dft(x::Vector{T}) where T <: Number
    N = length(x)

    X = zeros(Complex, N)

    for k in 0:N-1
        X[k+1] = sum([x[n+1] * exp(-2pi * im * k * n / N) for n in 0:N-1]) 
    end

    return X
end

function fft(x::Vector{T}) where T <: Number
    N = length(x)

    if ispow2(N)
        x_even = x[1:2:end]
        x_odd = x[2:2:end]

        factors = @. exp(-2pi * im * (0:N-1) / N)

        return vcat(x_even .+ x_odd .* factors[1:N÷2], x_even .+ x_odd .* factors[1+(N÷2):end])
    elseif isprime(N)
        return rader_fft(x)
    else
        factorisation = factor(N)
        factorisation[1]
    end
end

function fft(x::Vector{T}) where T <: Real
    N = length(x)

    if ispow2(N)
        x_even = fft(x[1:2:end])
        x_odd = fft(x[2:2:end])

        factors_first = @. exp(-2pi * im * (0:N÷2) / N)
        factors_rest = reverse(@. conj(factors_first[2:end-1]))

        factors = vcat(factors_first, factors_rest)

        return vcat(x_even .+ x_odd .* factors[1:N÷2], x_even .+ x_odd .* factors[1+(N÷2):end])
    elseif isprime(N)
        
    else
        
    end
end

display(dft([1, 2, 3]))