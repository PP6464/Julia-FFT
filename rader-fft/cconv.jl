# Shift first element
function shift(vector::Vector{T}) where T <: Number
    last = vector[end]

    result = zeros(T, length(vector))

    result[2:end] = vector[1:end - 1]
    result[1] = last

    result
end

display(shift([1, 2, 3, 4]))