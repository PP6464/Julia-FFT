include("../maths.jl")

function format_for_display(x::Union{Vector{T}, Matrix{T}}) where T <: Real
    return map(z -> round(z, digits = 2), x)
end

function format_for_display(x::Union{Vector{T}, Matrix{T}}) where T <: Number
    return map(z -> round(real(z), digits = 2) + im * round(imag(z), digits = 2), x)
end
function radix2FFT(x::Vector{T}, ω::Union{Complex{Float64}, Nothing} = nothing) where T <: Number
    N = length(x)

    if N == 1
        return x # One-point time domain spectrum is the same in frequency domain
    end

    #println("FFT of:")
    #display(x)
    #println("with root of unity: $ω")

    # Initialise the root of unity to be used (normally conjugate of principal Nᵗʰ root of unity)
    # (Can supply different ω values, see radix2IFFT)
    if ω === nothing
        ω = exp(-τ*im/N)
    end

    Xₑ = radix2FFT(x[1:2:end], ω^2) # Get the even indices of x (0-indexed) and recursively FFT (Julia is 1-indexed)
    Xₒ = radix2FFT(x[2:2:end], ω^2) # Get the odd indices of x (0-indexed) and recursively FFT (Julia is 1-indexed)

    res = zeros(Complex{Float64}, N)

    #println("X evens: ")
    #println(Xₑ)
    #println("X odds: ")
    #println(Xₒ)

    # Use the fact that the second part of the twiddle factors is the conjugate of the first when input is purely real
    for i in 1:N÷2
        #println("res[$i]: ")
        #println(Xₑ[i] + Xₒ[i] * ω^(i-1))
        #println("res[$(i+N÷2)]: ")
        #println(Xₑ[i] - Xₒ[i] * ω^(i-1))
        res[i] = Xₑ[i] + Xₒ[i] * ω^(i-1) 
        res[i+N÷2] = Xₑ[i] - Xₒ[i] * ω^(i-1)
    end

    return res
end

function radix2IFFT(x::Vector{T}, ω::Union{Complex{Float64}, Nothing} = nothing) where T <: Number
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
    u = [x[i+1] * ω^(-(i^2)/2) for i in 0:N-1]
    v = [ω^((i^2)/2) for i in 0:N-1]
    v_conj = [ω^(-(i^2)/2) for i in 0:N-1]

    L = 2^ceil(log2(2*N+1))
    ωₗ = exp(τ*im/L) # Calculate root of unity required

    uₗ = vcat(u, [0 for _ in 1:L-N]) # Pad u(n)

    # Pad v(n) and shift
    aux = v[2:end]
    reverse!(aux)
    vₗ = vcat(v, [0 for i in 1:L - 2*N + 1], aux)

    # FFT circular convolve
    fft_uₗ = radix2FFT(uₗ, ωₗ)
    fft_vₗ = radix2FFT(vₗ, ωₗ)

    uvₗ = fft_uₗ .* fft_vₗ

    iftₗ = radix2IFFT(uvₗ, ωₗ)

    ift = iftₗ[begin:N]

    # Multiply the IFT by v_conj to get result
    return v_conj .* ift
end

display(bluestein([i for i in 1:11]))
j = im
#display(radix2IFFT([(0.4999999999999987+0.8660254037844383j), (3.04088222241137-3.0355339059327404j), (11.232050807568882+4.133974596215561j), (5.408607520371811-4.0355339059327395j), (16.5+0.8660254037844384j), (-2.2370346451180003+4.035533905932738j), (7.767949192431123-5.866025403784442j), (5.787544902334822+3.0355339059327404j)], (0.7071067811865476+0.7071067811865476j)))