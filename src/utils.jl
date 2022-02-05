# Define custom function for nice imports
Base.parse(::Type{Vector{String}}, x::String) = Vector{String}(filter(x-> x != ", ", split(x, "\""))[2:end-1])
function Base.parse(::Type{Vector{Float64}}, x::String)
    number = split(split(x, "[")[end][1:end-1], ", ")
    number_list = Float64[]
    for num in number
        if num != ""
            push!(number_list, parse(Float64, num))
        else
            return push!(number_list, NaN)
        end
    end
    return number_list

end

function Base.parse(::Type{Tuple{Int64, Tuple{Int64, Int64}}}, x::String)
    numbers = split(x, ", ")
    ind1 = parse(Int64, split(numbers[1], "(")[2])
    ind2 = parse(Int64, split(numbers[2], "(")[2])
    ind3 = parse(Int64, split(numbers[3], ")")[1])
    return (ind1, (ind2, ind3))
end

Base.parse(::Type{Vector{String}}, x::Missing) = String[]
Base.parse(::Type{Vector{Float64}}, x::Missing) = Float64[]