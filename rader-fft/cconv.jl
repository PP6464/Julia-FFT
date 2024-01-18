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

function convert_to_matrix(v::Vector{Vector{T}}, columns::Bool = false) where T <: Number
    if columns
        m, n = length(v[1]), length(v)
        matrix = zeros(T, m, n)
        for i in 1:n
            matrix[begin:end, i] = v[i]
        end
        return matrix
    else 
        m, n = length(v), length(v[1])
        matrix = zeros(T, m, n)
        for i in 1:m
            matrix[i, begin:end] = v[i]
        end
        return matrix
    end
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
    if length(x) != length(y)
        throw(ArgumentError("The lists must be of the same length"))
    end

    N = length(x)

    primary_matrix = zeros(Int, N, N)

    for i in 1:N
        primary_matrix[1, i] = get(x, i, 0)
    end

    for i in 2:N
        primary_matrix[i, begin:end] = shift(rows(primary_matrix)[i - 1])
    end

    ultimate_matrix = convert_to_matrix([i for i in [transpose(primary_matrix)[i, begin:end] for i in 1:N]])

    return ultimate_matrix * y
end

# using BenchmarkTools

x = rand(Int, 20000)
y = rand(Int, 20000)

# x = [2 - i for i in 1:15]
# y = [i for i in 1:15]

@time display(cconv(x, y))