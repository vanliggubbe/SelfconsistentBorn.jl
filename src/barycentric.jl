abstract type RationalInterpolation <: Function end

mutable struct BarycentricInterpolation{X <: Number, W <: Number, T, M, A <: AbstractArray{T, M}} <: RationalInterpolation
    nodes :: Vector{X}
    weights :: Vector{W}
    values :: A
    perm :: Vector{Int}

    function BarycentricInterpolation(nodes, weights, values)
        if length(nodes) != length(weights)
            throw(DimensionMismatch("Arguments `nodes` and `weights` must have the same lengths"))
        end
        if length(nodes) > length(values)
            throw(DimensionMismatch("Array `values` must not be shorter then `nodes`"))
        end
        perm = collect(1 : length(weights))
        sortperm!(perm, abs.(weights))
        new{eltype(nodes), eltype(weights), eltype(values), ndims(values), typeof(values)}(
            nodes isa Vector ? nodes : vec(collect(nodes)),
            weights isa Vector ? weights : vec(collect(weights)),
            values,
            perm
        )
    end
end

function BarycentricInterpolation(nodes, weights, values :: AbstractVector{<: Array{T, N}}) where {T, N}
    if isempty(values)
        throw(ArgumentError("values must not be empty"))
    end
    sz = size(first(values))
    vals = Array{T, N + 1}(undef, (size(first(values))..., length(values)))
    
    for (i, val) in enumerate(values)
        if size(val) != sz
            throw(DimensionMismatch("All the values should have the same dimensions"))
        end
        selectdim(vals, N + 1, i) .= val
    end
    return BarycentricInterpolation(nodes, weights, vals)
end

evaluate!(_, :: BarycentricInterpolation{X, W, T, 1, A}, __) where {X, W, T, A} = error("Not applicable")
function evaluate!(y, F :: BarycentricInterpolation{X, W, T, M, A}, x) where {X, W, T, M, A}
    if size(y) != size(F.values)[begin : end - 1]
        throw(DimensionMismatch("Output array dimension mismatch: expected $(size(F.values)[begin : end - 1]), got $(size(y))"))
    end
    j = findfirst(let x = x; z -> z == x end, F.nodes)
    sz = size(y)
    len = prod(sz)
    if j isa Nothing
        fill!(y, zero(eltype(y)))
        # Float64 because of division
        den = zero(promote_type(X, W, typeof(x), Float64))
        tmp = zero(den)
        @inbounds for i in F.perm
            tmp = F.weights[i] / (x - F.nodes[i])
            den += tmp
            axpy!(tmp, selectdim(F.values, M, i), y)
        end
        return (y ./= den)
    else
        y .= selectdim(F.values, M, j)
        return y
    end
end

function (F :: BarycentricInterpolation{X, W, T, 1, A})(x) where {X, W, T, A}
    j = findfirst(let x = x; y -> y == x end, F.nodes)
    if j isa Nothing
        # Float64 because of division
        num = zero(promote_type(X, W, T, typeof(x), Float64))
        den = zero(promote_type(X, W, typeof(x), Float64))
        tmp = zero(den)
        @inbounds @simd for i in F.perm
            tmp = F.weights[i] / (x - F.nodes[i])
            num += tmp * F.values[i]
            den += tmp
        end
        return num / den
    else
        return F.values[j]
    end
end

function (F :: BarycentricInterpolation{X, W, T, N, A})(x) where {X, W, T, N, A}
    sz = size(F.values)[begin : end - 1]
    # Float64 because of division
    y = Array{promote_type(X, W, T, typeof(x), Float64)}(undef, sz)
    evaluate!(y, F, x)
end

function add_node!(F, node, weight, value)
    push!(F.nodes, node)
    push!(F.weights, weight)
    if length(F.nodes) > size(F.values, ndims(F.values))
        sz = size(F.values)
        tmp = vec(F.values)
        resize!(tmp, 2 * length(tmp) + prod(sz[begin : end - 1]))
        F.values = reshape(tmp, (sz[begin : end - 1]..., 2 * sz[end] + 1))
    end
    selectdim(F.values, ndims(F.values), length(F.nodes)) .= value
    return F
end

add_node!(F, node, value) = add_node!(F, node, zero(eltype(F.weights)), value)
