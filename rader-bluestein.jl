τ = 2π

using Primes

function format_for_display(x::Union{Vector{T}, Matrix{T}}) where T <: Real
    return map(z -> round(z, digits = 2), x)
end

function format_for_display(x::Union{Vector{T}, Matrix{T}}) where T <: Number
    return map(z -> round(real(z), digits = 2) + im * round(imag(z), digits = 2), x)
end

function radix2FFT(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    if N == 1
        return x # One-point time domain spectrum is the same in frequency domain
    end

    # Initialise the root of unity to be used (normally conjugate of principal Nᵗʰ root of unity)
    # (Can supply different ω values, see radix2IFFT)
    if ω === nothing
        ω = exp(-τ*im/N)
    end

    Xₑ = radix2FFT(x[1:2:end], ω^2) # Get the even indices of x (0-indexed) and recursively FFT (Julia is 1-indexed)
    Xₒ = radix2FFT(x[2:2:end], ω^2) # Get the odd indices of x (0-indexed) and recursively FFT (Julia is 1-indexed)

    res = zeros(ComplexF64, N)

    # Use the fact that the second part of the twiddle factors is the conjugate of the first when input is purely real
    for i ∈ 1:N÷2
        res[i] = Xₑ[i] + Xₒ[i] * ω^(i-1) 
        res[i+N÷2] = Xₑ[i] - Xₒ[i] * ω^(i-1)
    end

    return res
end

function radix2IFFT(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    if N == 1
        return x # One-point time domain spectrum is the same in frequency domain
    end

    # Initialise root of unity to be used (normally for IFFT this is principal Nᵗʰ root of unity) 
    if ω === nothing
        ω = exp(τ*im/N)
    end

    return radix2FFT(x, 1/ω) ./ N # Supply reciprocal of principal Nᵗʰ root of unity as IFFT = FFT but with that used instead to generate twiddle factors
end

function bluestein_FFT(x::Vector{T}, ω::Union{Complex{A}, Nothing} = nothing) where {T <: Number, A <: Number}
    N = length(x)

    if N == 1
        return x
    end

    # Initialise ω to exp(iτ/N) if not supplied
    # Supply different ω for IFFT
    if ω === nothing 
        ω = exp(τ*im/N)
    end

    # Make lists u(n) and v(n)
    u = [x[i+1] * ω^(-(i^2)/2) for i ∈ 0:N-1]
    v = [ω^((i^2)/2) for i ∈ 0:N-1]
    v_conj = [ω^(-(i^2)/2) for i ∈ 0:N-1]

    L = 2^ceil(log2(2*N+1))
    ωₗ = exp(log(ω) * N/L) # exp(τ*im/L) # Calculate root of unity required

    uₗ = vcat(u, [0 for _ ∈ 1:L-N]) # Pad u(n)

    # Pad v(n) and shift
    aux = v[2:end]
    reverse!(aux)
    vₗ = vcat(v, [0 for i ∈ 1:L - 2*N + 1], aux)

    # FFT circular convolve
    fft_uₗ = radix2FFT(uₗ, ωₗ)
    fft_vₗ = radix2FFT(vₗ, ωₗ)

    uvₗ = fft_uₗ .* fft_vₗ

    iftₗ = radix2IFFT(uvₗ, ωₗ)

    ift = iftₗ[begin:N]

    # Multiply the IFT by v_conj to get result
    return v_conj .* ift
end

function bluestein_IFFT(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    if N == 1
        return x
    end

    if ω === nothing
        ω = exp(τ*im/N)
    end

    return bluestein_FFT(x, 1/ω) ./ N
end

# Rader FFT

# Successive powers of g up to m - 1, all mod m
function sequence(g::Int64, m::Int64)
    return [(big(g)^i)%m for i ∈ 0:m-2]
end

function primroot(m::Integer)
    if !isprime(m)
        error("Argument 'm' must be a prime number")
    end
    if m == 2; return 1; end

    P = keys(factor(m-1))
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

# Use to generate the entire group (𝐙/p𝐙)* and successive powers (p is prime)
function generator(p::Int64)
    # for g ∈ 2:p-1
    #     perm = sequence(g, p)
    #     if length(perm) == p - 1
    #         return g, perm
    #     end
    # end
    return primroot(p), sequence(primroot(p), p)
end

# FFT convolution of u and v, where u and v are of the same length (uses padding)
function conv(u::Vector{T}, v::Vector{U}) where {T <: Number, U <: Number}
    N = length(u)
    L = 2^ceil(log2(2N+1))
    ωₗ = exp(τ*im/L)

    fft_uₗ = radix2FFT(vcat(u, [0 for _ ∈ 1:L-N]), ωₗ)
    fft_vₗ = radix2FFT(vcat(v, [0 for _ ∈ 1:L-N]), ωₗ)

    # Convolution theorem: 𝓕{u(x) ∗ v(x)} = 𝓕{u(x)} . 𝓕{v(x)}
    uₗvₗ = radix2IFFT(fft_uₗ .* fft_vₗ)

    return [uₗvₗ[i] + uₗvₗ[i + N] for i in 1:N]
end

function rader_FFT(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    # Initialise ω to exp(iτ/N) if not supplied
    # See rader_IFFT for when different ω supplied
    if ω === nothing
        ω = exp(-τ*im/N)
    end

    g, g_seq = generator(N)
    g⁻¹ = (g^(N-2))%N 
    g⁻¹_seq = sequence(g⁻¹, N)

    # u(n) = ω^g⁻ⁿ
    u = [ω^gₙ for gₙ ∈ g⁻¹_seq]

    # v(n) = x[gᵐ]
    v = [x[gᵐ+1] for gᵐ ∈ g_seq]

    uv = conv(u, v)

    fft = zeros(ComplexF64, N)
    fft[1] = ∑(x) # FFT[0] = ∑ᵢ xᵢ (Julia is 1-indexed)

    # FFT[g⁻ʲ] = a₀ + (u * v)ⱼ (Julia is 1-indexed)
    for i in 1:N-1
        fft[g⁻¹_seq[i] + 1] = x[1] + uv[i]
    end

    return fft
end

function rader_IFFT(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    # Initialise ω to principal Nth root of unity if not supplied
    if ω === nothing
        ω = exp(τ*im/N)
    end

    return rader_FFT(x, 1/ω) ./ N
end

function ispow2(x::Int64)
    return (x & x-1) == 0 # Bitwise and of x and x-1 == 0 ⟺ x is a power of 2
end

function fft(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    # Initialise ω to exp(iτ/N) if not supplied
    # See IFFT for when different ω supplied
    if ω === nothing
        ω = exp(im*τ/N)
    end

    if ispow2(N)
        return radix2FFT(x, ω)
    else
        return bluestein_FFT(x, ω)
    end
end

function ifft(x::Vector{T}, ω::Union{Complex{U}, Nothing} = nothing) where {T <: Number, U <: Number}
    N = length(x)

    # Initialise ω to principal Nth root of unity if not supplied
    if ω === nothing
        ω = exp(τ*im/N)
    end

    return fft(x, 1/ω) ./ N
end

function fft_convolve(x::Vector{T}, y::Vector{U}) where {T <: Number, U <: Number}
    X = vcat(x, zeros(T, length(y) - 1))
    Y = vcat(y, zeros(U, length(x) - 1))

    return ifft(fft(X) .* fft(Y))[1:length(x) + length(y) - 1]
end

# Convert Matrix{T} into Vector{Vector{T}} (a vector of the columns of the input matrix)
function columns(x::Matrix{T})::Vector{Vector{T}} where T <: Any
    return [x[begin:end, i] for i ∈ 1:size(x)[2]]
end

# Convert Matrix{T} into Vector{Vector{T}} (a vector of the rows of the input matrix)
function rows(x::Matrix{T})::Vector{Vector{T}} where T <: Any
    return [x[i, begin:end] for i ∈ 1:size(x)[1]]
end

# Converts Vector{Vector{T}} into Matrix{T}. Columns specifies whether the input is a vector of columns (true) or a vector of rows (false)
function convert_to_matrix(x::Vector{Vector{T}}, columns::Bool = false)::Matrix{T} where T <: Any
    if columns
        res = zeros(T, length(x[1]), length(x))
        for i ∈ 1:length(x)
            res[begin:end, i] = x[i]
        end
        return res
    else
        res = zeros(T, length(x), length(x[1]))
        for i ∈ 1:length(x)
            res[i, begin:end] = x[i]
        end
        return res
    end
end

function fft2(x::Matrix{T}, ω::Union{Complex{A}, Nothing} = nothing)::Matrix{ComplexF64} where {T <: Number, A <: Number}
    fft_columns = [fft(i, ω) for i ∈ columns(x)]
    fft_rows = [fft(i, ω) for i ∈ rows(convert_to_matrix(fft_columns, true))]
    return convert_to_matrix(fft_rows)
end

function ifft2(x::Matrix{T}, ω::Union{Complex{A}, Nothing} = nothing)::Matrix{ComplexF64} where {T <: Number, A <: Number}
    ifft_columns = [ifft(i, ω) for i ∈ columns(x)]
    ifft_rows = [ifft(i, ω) for i ∈ rows(convert_to_matrix(ifft_columns, true))]
    return convert_to_matrix(ifft_rows)
end

function conv2(x::Matrix{T}, y::Matrix{U})::Matrix{ComplexF64} where {T <: Number, U <: Number}
    # Pad to correct dimensions
    new_m = size(x)[1] + size(y)[1] - 1
    new_n = size(x)[2] + size(y)[2] - 1

    new_x = zeros(T, new_m, new_n)
    new_y = zeros(U, new_m, new_n)
    
    new_x[1:size(x)[1], 1:size(x)[2]] = x
    new_y[1:size(y)[1], 1:size(y)[2]] = y

    x = new_x
    y = new_y
    
    return ifft2(fft2(x) .* fft2(y))
end

# l = 2^20
# using BenchmarkTools

# x = [rand() for _ in 1:l]

# @time begin
#     display("RADIX 2: ")
#     radix2FFT(x)
# end

# @time begin
#     display("BLUESTEIN: ")
#     bluestein_FFT(x)
# end

function oaconvolve(x::Vector{Int64}, y::Vector{Int64}, block_size::Int64 = 4)
    # Are x and y's lengths mulitples of block_size?
    X = length(x)
    Y = length(y)

    # In case the lists are padded
    orig_X = X
    orig_Y = Y

    x_mul_block_size = false
    y_mul_block_size = false

    if X % block_size == 0
        x_mul_block_size = true
    elseif Y % block_size == 0
        y_mul_block_size = true
    end

    chunks = [[] for _ ∈ 1:(X>Y ? X : Y)]
    divided_signal = "x"

    if x_mul_block_size
        # Divide x into chunks
        for i ∈ 1:(X÷block_size)
            chunks[i] = x[(i-1) * block_size + 1:i * block_size]
        end
    elseif y_mul_block_size
        # Divide y into chunks
        divided_signal = "y"
        for i ∈ 1:(Y÷block_size)
            chunks[i] = y[(i-1) * block_size + 1:i * block_size]
        end
    else
        # Pad x and then divide it into pieces
        n_zeros = block_size - (X % block_size)
        x = vcat(x, zeros(Int64, n_zeros))
        X = block_size * ((X÷block_size) + 1) # Reflect change in x
        for i ∈ 1:(X÷block_size)
            chunks[i] = x[(i-1) * block_size + 1:i * block_size]
        end
    end

    while !isempty(chunks) && chunks[end] == []
        pop!(chunks)
    end

    output_chunks = [[] for _ ∈ eachindex(chunks)]

    for i ∈ enumerate(chunks)
        if divided_signal == "x"
            output_chunks[i[1]] = fft_convolve(convert(Vector{Int64}, map(j -> round(real(j)), i[2])), y)
        else
            output_chunks[i[1]] = fft_convolve(convert(Vector{Int64}, map(j -> round(real(j)), i[2])), x)
        end
    end

    desired_length = X+Y-1

    # Shift the output chunks to make reconstruction easier
    for i ∈ eachindex(chunks)
        if (desired_length - length(output_chunks[i]) - (i-1)*block_size) < 0
            println(desired_length)
            println(length(output_chunks[i]))
            println(block_size)
            println((i - 1) * block_size)
        end
        output_chunks[i] = vcat(
            zeros(Int64, (i-1) * block_size),
            output_chunks[i],
            zeros(Int64, desired_length - length(output_chunks[i]) - (i-1)*block_size)
        )
    end

    output = zeros(BigInt, X+Y-1)
    for i ∈ 1:desired_length
        output[i] = sum(map(x -> round(real(x[i])), output_chunks))
    end

    return output[1:(orig_X + orig_Y - 1)]
end

function oaconvolve2(x::Matrix{T}, y::Matrix{U}, block_size_row::Int64 = 4, block_size_col::Int64 = 4) where {T, U <: Number}
    x_rows, x_cols, y_rows, y_cols = size(x)[1], size(x)[2], size(y)[1], size(y)[2]
    orig_x_rows, orig_x_cols, orig_y_rows, orig_y_cols = x_rows, x_cols, y_rows, y_cols # In case a matrix is padded

    # Decide whether to divide up x or y
    divided_mat = "x"
    chunks = [] # Chunks stored as a list of rows

    if x_rows % block_size_row == 0 && x_cols % block_size_col == 0
        # x shall be divided, so let us chunk it
        chunks = [x[(i-1)*block_size_row+1:i*block_size_row, begin:end] for i ∈ 1:(x_rows÷block_size_row)]
        chunks = [[chunks[i][begin:end, (j-1)*block_size_col+1:j*block_size_col] for j ∈ (1:x_cols÷block_size_col)] for i ∈ 1:(x_rows÷block_size_row)]
        chunks = reduce(vcat, chunks) # Flatten to a vector of matrices
    elseif y_rows % block_size_row == 0 && y_cols % block_size_col == 0
        # y shall be divided, so let us chunk it
        divided_mat = "y"
        chunks = [y[(i-1)*block_size_row+1:i*block_size_row, begin:end] for i ∈ 1:(y_rows÷block_size_row)]
        chunks = [[chunks[i][begin:end, (j-1)*block_size_col+1:j*block_size_col] for j ∈ (1:y_cols÷block_size_col)] for i ∈ 1:(x_rows÷block_size_row)]
        chunks = reduce(vcat, chunks) # Flatten to a vector of matrices
    else
        # x will be divided, but needs to be padded to a nice size
        padded_x = zeros(T, x_rows % block_size_row == 0 ? x_rows : x_rows + (block_size_row - x_rows % block_size_row), x_cols % block_size_col == 0 ? x_cols : x_cols + (block_size_col - x_cols % block_size_col))
        padded_x[1:x_rows, 1:x_cols] = x
        x = padded_x
        x_rows, x_cols = size(x)
        chunks = [x[(i-1)*block_size_row+1:i*block_size_row, begin:end] for i ∈ 1:(x_rows÷block_size_row)]
        chunks = [[chunks[i][begin:end, (j-1)*block_size_col+1:j*block_size_col] for j ∈ (1:x_cols÷block_size_col)] for i ∈ 1:(x_rows÷block_size_row)]
        chunks = reduce(vcat, chunks) # Flatten to a vector of matrices
    end

    # Convolve each chunk with the other matrix
    chunks = map(m -> conv2(m, divided_mat == "x" ? y : x), chunks)

    # Shift each chunk into the correct position
    chunks = map(enumerate(chunks)) do i
        index = i[1]
        chunk = i[2]
        shifted_mat = zeros(Number, x_rows + y_rows - 1, x_cols + y_cols - 1)
        column_shifted_to = index % ((divided_mat == "x" ? x_cols : y_cols) ÷ block_size_col)
        row_shifted_to = index ÷ ((divided_mat == "x" ? x_cols : y_cols) ÷ block_size_col)
        shifted_mat[row_shifted_to*block_size_row+1:row_shifted_to*block_size_row+size(chunk)[1], column_shifted_to*block_size_col+1:column_shifted_to*block_size_col+size(chunk)[2]] = chunk
        return shifted_mat
    end

    output = zeros(Number, x_rows + y_rows - 1, x_cols + y_cols - 1)

    for i ∈ 1:x_cols+y_cols-1
        for j ∈ 1:x_rows+y_rows-1
            output[i,j] = sum(map(c -> c[i,j], chunks))
        end
    end

    return output[1:orig_x_rows+orig_y_rows-1, 1:orig_x_cols+orig_y_cols-1]
end

# l = 2*20

# list1 = rand(Int, l)
# list2 = rand(Int, l)

# using BenchmarkTools

# @time begin
# println("FFT convolve")
# fft_convolve(list1, list2)
# end

# @time begin
# println("OA convolve")
# oaconvolve(list1, list2, 16)
# end

oaconvolve2([1 2 3 4;5 6 7 8;9 10 11 12;13 14 15 16], [1 2 3 4;3 4 5 6;5 6 7 8;2 3 4 5], 2, 3)