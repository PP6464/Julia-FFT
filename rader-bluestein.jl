include("../maths.jl")

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

function bluestein(x::Vector{T}) where T <: Number
    N = length(x)

    ω = exp(τ*im/N)

    # Make lists u(n) and v(n)
    u = [x[i+1] * ω^(-(i^2)/2) for i ∈ 0:N-1]
    v = [ω^((i^2)/2) for i ∈ 0:N-1]
    v_conj = [ω^(-(i^2)/2) for i ∈ 0:N-1]

    L = 2^ceil(log2(2*N+1))
    ωₗ = exp(τ*im/L) # Calculate root of unity required

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

# FFT convolution of u and v, where u and V are of the same length (uses padding)
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

function rader(x::Vector{T}) where T <: Number
    N = length(x)

    g, g_seq = generator(N)
    g⁻¹ = (g^(N-2))%N 
    g⁻¹_seq = sequence(g⁻¹, N)

    # u(n) = ω^g⁻ⁿ
    ωₚ = exp(τ*im/N) # N-th root of unity
    u = [ωₚ^gₙ for gₙ ∈ g⁻¹_seq]

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

j = im

display(rader([1im,2,3]))