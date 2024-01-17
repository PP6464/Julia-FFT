using FFTW

function nextpow2(x::Int64)
    return 2 ^ ceil(log2(x))
end

function bluestein_fft(x::Vector)
    n = length(x)

    m = nextpow2(2n-1)
    y = vcat(x, zeros(ComplexF64, m - n))

    k = 0:m-1
    ω = exp.(-2pi * im /(2m) * k .^ 2)
end