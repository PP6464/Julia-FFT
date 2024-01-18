include("../maths.jl")

function radix2FFT(x::Vector{T}) where T <: Number
    N = length(x)

    if ispow2(N)
        x_even = x[1:2:end]
        x_odd = x[2:2:end]

        factors = @. exp(-2pi * im * (0:N-1) / N)

        return vcat(x_even .+ x_odd .* factors[1:N÷2], x_even .+ x_odd .* factors[1+(N÷2):end])
    else
        throw(ArgumentError("Length must be a power of 2"))
    end
end

function radix2IFFT(x::Vector{T}, normalise::Bool = true) where T <: Number
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
end 