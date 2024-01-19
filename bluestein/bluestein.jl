include("../maths.jl")

function radix2FFT(x::Vector{Int64})
    N = length(x)

    if !ispow2(N)
        throw(ArgumentError("List must have a length that is a power of 2"))
    end

end