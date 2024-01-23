include("maths.jl")

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

    res = zeros(Complex{Float64}, N)

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
    return [(g^i)%m for i ∈ 0:m-2]
end

# Use to generate the entire group (𝐙/p𝐙)* and successive powers (p is prime)
function generator(p::Int64)
    for g ∈ 2:p-1
        perm = sequence(g, p)
        if length(perm) == p - 1
            return g, perm
        end
    end
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
        ω = exp(τ*im/N)
    end

    g, g_seq = generator(N)
    g⁻¹ = (g^(N-2))%N 
    g⁻¹_seq = sequence(g⁻¹, N)

    # u(n) = ω^g⁻ⁿ
    u = [ω^gₙ for gₙ ∈ g⁻¹_seq]

    # v(n) = x[gᵐ]
    v = [x[gᵐ+1] for gᵐ ∈ g_seq]

    uv = conv(u, v)

    fft = zeros(Complex{Float64}, N)
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

function int_to_bits(x::Int64, total_length::Int64)
    return digits(x, base = 2, pad = total_length)
end

function bits_to_int(x::Vector{Int64})
    return sum([v*2^(i-1) for (i, v) ∈ enumerate(x)])
end

# Computes the FFT butterfly operator to reconstruct an 2N-point spectrum from two N point spectra
function fft_butterfly(evens::Vector{A}, odds::Vector{B}, ω::Union{Complex{C}, Nothing} = nothing) where {A <: Number, B <: Number, C <: Number}
    N = length(evens)

    # Initialise ω to exp(-iτ/2N)
    if ω === nothing
        ω = exp(im*τ/2N)
    end

    res = zeros(Complex, 2N)

    for i ∈ 1:N
        res[i] = evens[i] + odds[i] * ω^(i-1)
        res[i+N] = evens[i] - odds[i] * ω^(i-1)
    end

    return res
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
    elseif isprime(N)
        return rader_FFT(x, ω)
    else
        factors = factor(N)

        if 2 ∈ keys(factors) && length(factors) == 2 && factors[collect(keys(factors))[end]] == 1 # N can be expressed in the form 2ᵏ ⋅ p, where k ∈ ℕ, p is prime
            k = factors[2]
            p = collect(keys(factors))[end]

            lₚ = length(digits(p, base=2)) # Length of p in base 2

            # Can split into 2ᵏ chunks of size p, then rader FFT each chunk, then reconstruct
            chunks = collect(Iterators.partition(x, p))
            chunk_indices = [(i-1) for i ∈ 1:p]
            chunk_bit_indices = [int_to_bits(i, lₚ) for i ∈ chunk_indices]
            chunk_bit_reversed_indices = [bits_to_int(reverse(i)) + 1 for i ∈ chunk_bit_indices] # Julia is 1-indexed, hence adding 1 after bit reversal
            chunks_rearranged = chunks[chunk_bit_reversed_indices]
            chunks_fft = [rader_FFT(chunk) for chunk ∈ chunks_rearranged]

            res::Vector{Vector{Complex{Float64}}} = chunks_fft

            for i ∈ 1:k
                butterfly_size = 2^(i-1)
                res_chunked = collect(Iterators.partition(res_copy, butterfly_size))
                res = [fft_butterfly(i[1], i[2], ω) for i ∈ res_chunked]
            end

            return res
        else

        end
    end
end

display(fft([1,2,3,4,5,6]))