using LinearAlgebra

# Returns a vector of the rows of the matrix
function rows(matrix::Matrix)
    m = size(matrix)[begin]

    return [matrix[i, begin:end] for i in 1:m]
end

# Returns a vector of the columns of the matrix
function columns(matrix::Matrix)
    n = size(matrix)[end]

    return [matrix[begin:end, i] for i in 1:n]
end

# Shift first element
function shift(vector::Vector{T}) where T <: Number
    last = vector[end]

    result = zeros(T, length(vector))

    result[2:end] = vector[1:end - 1]
    result[1] = last

    result
end

function cconv(x::Vector{Int64}, y::Vector{Int64})
    N, M = length(x), length(y)

    size_of_matrix = max(N, M)

    primary_matrix = zeros(Int, size_of_matrix, size_of_matrix)

    for i in 1:size_of_matrix
        primary_matrix[1, i] = x[i]
    end

    for i in 2:size_of_matrix
        primary_matrix[i, begin:end] = shift(rows(primary_matrix)[i - 1])
    end

    display(primary_matrix)

    ultimate_matrix = convert_to_matrix([i for i in rows(transpose(primary_matrix))])
    difference_in_length = abs(N - M)

    for _ in M:M+difference_in_length - 1
        push!(y, 0)
    end

    return dot(ultimate_matrix, y)
end

display(cconv([2, 1, 2, 1], [1, 2, 3, 4]))