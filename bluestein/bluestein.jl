include("../maths.jl")

function radix2FFT(x::Vector{T}, ω::Union{Complex{Float64}, Nothing} = nothing) where T <: Number
    N = length(x)

    if N == 1
        return x
    end

    if ω === nothing
        ω = exp(-τ*im/N)
    end

    ω² = ω^2
    Xₑ = radix2FFT(x[1:2:end], ω²)
    Xₒ = radix2FFT(x[2:2:end], ω²)

    res = zeros(Complex{Float64}, N)

    for i in 1:N÷2
        res[i] = Xₑ[i] + Xₒ[i] * ω^i 
        res[i+N÷2] = Xₑ[i] - Xₒ[i] * ω^i
    end

    return res
end

function radix2IFFT(x::Vector{T}, ω::Union{Complex{Float64}, Nothing} = nothing) where T <: Number
    N = length(x)

    if ω === nothing
        ω = exp(τ*im/N)
    end

    if N == 1
        return x
    end

    return radix2FFT(x, 1/ω) ./ N

end

display(radix2FFT([1]))