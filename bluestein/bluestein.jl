include("../maths.jl")

function radix2FFT(x::Vector{T}) where T <: Number 
    N = length(x)

    if !ispow2(N)
        throw(ArgumentError("Ensure list is a length of a power of 2"))
    elseif N == 1
        return x
    end
    
    X_even = radix2FFT(x[1:2:end])
    X_odd = radix2FFT(x[2:2:end])
    factors::Vector{ComplexF64} = [] 
    
    if T <: Real
        factors_first = @. exp(-2pi * im * (0:N÷2) / N)
        factors_rest = reverse(@. conj(factors_first[2:end - 1]))
        factors = vcat(factors_first, factors_rest)
    else
        factors = @. exp(-2pi * im * (0:N-1) / N)
    end

    return vcat(X_even .+ factors[begin:N÷2] .* X_odd, X_even .+ factors[(1 + N÷2):end] .* X_odd)
end

function ifft(x::Vector{T}, normalise::Bool = true) where T <: Number
    N = length(x)

    if !ispow2(N)
        throw(ArugmentError("Ensure list is a length of a power of 2"))
    elseif N == 1
        return x
    end

    X_even = ifft(x[1:2:end])
    X_odd = ifft(x[2:2:end])
    factors = @. exp(2pi * im * (0:N-1) / N)

    return vcat(X_even .+ factors[begin:N÷2] .* X_odd, X_even .+ factors[1 + N÷2:end] .* X_odd)  ./ (normalise ? sqrt(N) : 1)
end

function bluestein_FFT(x::Vector{T}) where T <: Number
    N = length(x)
    ω = exp(-2pi*im*𝜏)

    u = [x[i] * ω^(-(i^2)/2) for i in 1:N]
    v = [ω^((i^2)/2) for i in 1:N]
    v_conj = @. conj(v)

    L = 2^ceil(log2(2N + 1))
    ω_l = exp(𝜏/L)

    u_padded = vcat(u, [0 for _ in 1:L])

    aux = v[2:end]
    reverse!(aux)
    v_padded = vcat(v, [0 for _ in 1:(L - 2N + 1)], aux)

    fft_u_padded = radix2FFT(u_padded)
    fft_v_padded = radix2FFT(v_padded)

    uv_padded = fft_u_padded .* fft_v_padded

    ift_padded = radix2IFFT(uv_padded, ω_l)

    return ift_padded .* v_conj
end

display(bluestein_FFT([1, 2, 3, 4, 5]))