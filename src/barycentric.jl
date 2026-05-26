mutable struct BarycentricInterpolation{M, V, A <: ElasticArray{V, M}, X <: Number, W <: Number} <: RationalInterpolation{M}
    nodes :: Vector{X}
    weights :: Vector{W}
    values :: A
    perm :: Vector{Int}

    function BarycentricInterpolation(nodes, weights, values)
        if length(nodes) != length(weights) || size(values)[end] != length(nodes)
            throw(DimensionMismatch("Arguments `nodes` and `weights` must have the same lengths"))
        end
        perm = collect(eachindex(weights))
        sortperm!(perm, abs.(weights))
        new{ndims(values), eltype(values), typeof(values), eltype(nodes), eltype(weights)}(
            nodes isa Vector ? nodes : vec(collect(nodes)),
            weights isa Vector ? weights : vec(collect(weights)),
            values,
            perm
        )
    end
end

Base.ndims(:: BarycentricInterpolation{M}) where {M} = M
Base.ndims(:: Type{<: BarycentricInterpolation{M}}) where {M} = M

function BarycentricInterpolation(nodes, weights, values :: AbstractVector{<: Array{T, N}}) where {T, N}
    if isempty(values)
        throw(ArgumentError("values must not be empty"))
    end
    vals = reduce(
        append!,
        values;
        init = ElasticArray{T, N + 1}(undef, size(first(values))..., 0)
    )
    return BarycentricInterpolation(nodes, weights, vals)
end

evaluate!(_, :: BarycentricInterpolation{1}, __) = error("Not applicable")

function evaluate!(y, F :: BarycentricInterpolation, x :: Number)
    if size(y) != size(F.values)[begin : end - 1]
        throw(DimensionMismatch("Output array dimension mismatch: expected $(size(F.values)[begin : end - 1]), got $(size(y))"))
    end
    j = findfirst(let x = x; z -> isapprox(z, x) end, F.nodes)
    sz = size(y)
    if j isa Nothing
        fill!(y, zero(eltype(y)))
        den = zero(divtype(eltype(F.weights), promote_type(typeof(x), eltype(F.nodes))))
        tmp = zero(den)
        @inbounds for i in F.perm
            tmp = F.weights[i] / (x - F.nodes[i])
            den += tmp
            axpy!(tmp, slice(F.values, i), y)
        end
        return (y ./= den)
    else
        y .= slice(F.values, j)
        return y
    end
end

function (F :: BarycentricInterpolation{1})(x)
    j = findfirst(let x = x; y -> isapprox(x, y) end, F.nodes)
    if j isa Nothing
        num = zero(
            divtype(
                promote_type(eltype(F.weights), eltype(F.values)),
                promote_type(typeof(x), eltype(F.nodes))
            )
        )
        den = zero(divtype(eltype(F.weights), promote_type(typeof(x), eltype(F.nodes))))
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

function (F :: BarycentricInterpolation)(x)
    sz = size(F.values)[begin : end - 1]
    T = divtype(
        divtype(promote_type(eltype(F.values), eltype(F.weights)), promote_type(eltype(F.nodes), typeof(x))),
        divtype(eltype(F.weights), promote_type(eltype(F.nodes), typeof(x)))
    )
    y = Array{T}(undef, sz)
    evaluate!(y, F, x)
end

function add_node!(F, node, weight, value)
    if any(isapprox.(node, F.nodes))
        throw(ArgumentError("Node $(node) already exists"))
    end
    push!(F.nodes, node)
    push!(F.weights, weight)
    append!(F.values, value)
    push!(perm, lastindex(F.weights))
    sortperm!(perm, abs.(weights))
    return F
end

add_node!(F, node, value) = add_node!(F, node, zero(eltype(F.weights)), value)

mutable struct OddBarycentricInterpolation{M, V, A <: ElasticArray{V, M}, X <: Real, W <: Number} <: RationalInterpolation{M}
    nodes :: Vector{X}
    weights :: Vector{W}
    values :: A
    perm :: Vector{Int}

    function OddBarycentricInterpolation(nodes, weights, values)
        if length(nodes) != length(weights) || size(values)[end] != length(nodes)
            throw(DimensionMismatch("Arguments `nodes` and `weights` must have the same lengths"))
        end
        if any(nodes .<= 0.0)
            throw(ArgumentError("All the nodes must be positive"))
        end
        perm = collect(eachindex(weights))
        sortperm!(perm, abs.(weights))
        values′ = values isa ElasticArray ? values : ElasticArray(values)
        new{ndims(values), eltype(values), typeof(values′), eltype(nodes), eltype(weights)}(
            nodes isa Vector ? nodes : vec(collect(nodes)),
            weights isa Vector ? weights : vec(collect(weights)),
            values′,
            perm
        )
    end
end

function OddBarycentricInterpolation(nodes, weights, values :: AbstractVector{<: Array{T, N}}) where {T, N}
    if isempty(values)
        throw(ArgumentError("values must not be empty"))
    end
    vals = reduce(
        append!, values;
        init = ElasticArray{T}(undef, size(first(values))..., 0)
    )
    return OddBarycentricInterpolation(nodes, weights, vals)
end


# TODO implement a better version
(F :: OddBarycentricInterpolation)(x) = sum(
    slice(F.values, i) * F.weights[i] * F.nodes[i] / (x ^ 2 - F.nodes[i] ^ 2) for i in F.perm
) / sum(F.weights[i] * x / (x ^ 2 - F.nodes[i] ^ 2) for i in F.perm)

(F :: OddBarycentricInterpolation{1})(x) = sum(
    F.values[i] * F.weights[i] * F.nodes[i] / (x ^ 2 - F.nodes[i] ^ 2) for i in F.perm
) / sum(F.weights[i] * x / (x ^ 2 - F.nodes[i] ^ 2) for i in F.perm)

struct SymmetricBarycentricInterpolation{M, V, A <: ElasticArray{V, M}, X <: Number, W <: Number, SN, SW, SV} <: RationalInterpolation{M}
    nodes :: Vector{X}
    weights :: Vector{W}
    values :: A

    sym_n :: SN
    sym_w :: SW
    sym_v :: SV
    perm :: Vector{Int}

    function SymmetricBarycentricInterpolation(nodes, weights, values, sym_n, sym_w, sym_v)
        if length(weights) ≠ last(size(values))
            throw(DimensionMismatch("Number of weights must coincide with the number of values"))
        end
        if length(nodes) ∉ [length(weights), length(weights) - 1]
            throw(DimensionMismatch("Number of nodes must coincide with the number of weights or be less by one"))
        end
        if any(nodes .<= 0.0)
            throw(ArgumentError("All the nodes must be positive"))
        end
        perm = collect(eachindex(nodes))
        sortperm!(perm, abs.(weights[1 : length(nodes)]))
        values′ = values isa ElasticArray ? values : ElasticArray(values)
        new{ndims(values), eltype(values), typeof(values′), eltype(nodes), eltype(weights), typeof(sym_n), typeof(sym_w), typeof(sym_v)}(
            nodes isa Vector ? nodes : vec(collect(nodes)),
            weights isa Vector ? weights : vec(collect(weights)),
            values′,
            sym_n,
            sym_w,
            sym_v,
            perm
        )
    end
end

evaluate!(__, :: SymmetricBarycentricInterpolation{1}, _) = error("Not applicable")

function evaluate!(y, F :: SymmetricBarycentricInterpolation, x :: Number)
    if size(y) != firstdims(F.values)
        throw(DimensionMismatch("Output array dimension mismatch: expected $(size(F.values)[begin : end - 1]), got $(size(y))"))
    end
    j = findfirst(let x = x; z -> isapprox(z, x) end, F.nodes)
    if j isa Nothing
        j = findfirst(let x = -x; z -> isapprox(z, x) end, F.nodes)
        if j isa Nothing
            fill!(y, zero(eltype(y)))
            den = zero(divtype(eltype(F.weights), promote_type(typeof(x), eltype(F.nodes))))
            if length(F.nodes) < length(F.weights)
                y .= slice(F.values, size(F.values, ndims(F.values)))
                den = last(F.weights)
            end
            tmp1 = zero(den)
            tmp2 = zero(den)
            @inbounds for i in F.perm
                tmp1 = F.weights[i] / (x - F.nodes[i])
                tmp2 = F.sym_w(F.weights[i]) / (x - F.sym_n(F.nodes[i]))
                den += tmp1
                den += tmp2
                axpy!(tmp1, slice(F.values, i), y)
                evaluate!(y, F.sym_v, slice(F.values, i), tmp2, one(eltype(y)))
            end
            return (y ./= den)
        else
            evaluate!(y, F.sym_v, slice(F.values, j))
            return y
        end
    else
        y .= slice(F.values, j)
        return y
    end
end

function evaluate_diff!(y, F :: SymmetricBarycentricInterpolation, x)
    if size(y) != firstdims(F.values)
        throw(DimensionMismatch("Output array dimension mismatch: expected $(size(F.values)[begin : end - 1]), got $(size(y))"))
    end
    j = findfirst(let x = x; z -> isapprox(z, x) end, F.nodes)
    if j isa Nothing
        j = findfirst(let x = -x; z -> isapprox(z, x) end, F.nodes)
        if j isa Nothing
            fill!(y, zero(eltype(y)))
            den = zero(divtype(eltype(F.weights), promote_type(typeof(x), eltype(F.nodes))))
            tmp1 = zero(den)
            tmp2 = zero(den)
            @inbounds for i in F.perm
                tmp1 = F.weights[i] / (x - F.nodes[i])
                tmp2 = F.sym_w(F.weights[i]) / (x + F.nodes[i])
                den += tmp1
                den += tmp2
                axpy!(tmp1, slice(F.values, i), y)
                evaluate!(y, F.sym_v, slice(F.values, i), tmp2, one(eltype(y)))
            end
            return (y ./= den)
        else
            evaluate!(y, F.sym_v, slice(F.values, j))
        end
    else
        y .= slice(F.values, j)
        return y
    end
end


function (F :: SymmetricBarycentricInterpolation)(x)
    sz = size(F.values)[begin : end - 1]
    T = divtype(
        divtype(promote_type(eltype(F.values), eltype(F.weights)), promote_type(eltype(F.nodes), typeof(x))),
        divtype(eltype(F.weights), promote_type(eltype(F.nodes), typeof(x)))
    )
    y = Array{T}(undef, sz)
    evaluate!(y, F, x)
end

#=

sum_j w_j / (x - x_j) + w̄_j / (x + x_j) =
sum_j [w_j (x + x_j) + w̄_j (x - x_j)] / (x² - x_j²) = 
sum_j (2 Re w_j x + 2i Im w_j x_j) / (x² - x_j²)
=#

function poles(F :: SymmetricBarycentricInterpolation)
    B = Matrix(1.0I, 1 + 2 * length(F.nodes), 1 + 2 * length(F.nodes))
    B[end, end] = 0.0
    A = zeros(promote_type(real(eltype(F.weights)), eltype(F.nodes)), size(B))
    @inbounds for i in eachindex(F.weights, F.nodes)
        A[2 * i - 1, 2 * i - 0] = F.nodes[i]
        A[2 * i - 0, 2 * i - 1] = -F.nodes[i]
        A[end, 2 * i - 1] = 1.0
        A[2 * i - 1, end] = real(F.weights[i])
        A[2 * i - 0, end] = imag(F.weights[i])
        #A[1, 1 + i] = sqrt(abs(F.weights[i]))
        #A[1, 1 + i + length(F.weights)] = sqrt(abs(F.sym_w.(F.weights[i])))
        #A[1 + i, 1] = F.weights[i] / A[1, 1 + i]
        #A[1 + i + length(F.weights), 1] = F.sym_w.(F.weights[i]) / A[1, 1 + i + length(F.weights)]
        #A[1 + i, 1 + i] = F.nodes[i]
        #A[1 + i + length(F.weights), 1 + i + length(F.weights)] = -F.nodes[i]
    end
    return im * filter(isfinite, eigvals(A, B; sortby = real))
end

#=
function Base.write(parent :: Union{File, Group}, name :: AbstractString, f :: BarycentricInterpolation)
    g = create_group(parent, name)
    g["nodes"] = f.nodes
    g["weights"] = f.weights
    g["values"] = f.values

    attributes(g)["__julia_type__"] = String(typeof(f))
    return g
end

function Base.read(parent::Union{File, Group}, name :: AbstractString, ::Type{<: BarycentricInterpolation})
    g = parent[name]
    return BarycentricInterpolation(
        read(g, "nodes"),
        read(g, "weights"),
        ElasticArray(read(g, "values"))
    )
end

function Base.read(parent::Union{File, Group}, name :: AbstractString, :: Type{<: SymmetricBarycentricInterpolation})
    g = parent[name]
    return SymmetricBarycentricInterpolation(
        read(g, "nodes"),
        read(g, "weights"),
        ElasticArray(read(g, "values")),
        conj,
        conj ∘ PermutationMap(read(g, "sym_perm"))
    )
end

function Base.write(parent :: Union{File, Group}, name :: AbstractString, f :: SymmetricBarycentricInterpolation)
    g = create_group(parent, name)
    g["nodes"] = f.nodes
    g["weights"] = f.weights
    g["values"] = f.values
    g["sym_perm"] = f.sym_v.inner.p

    attributes(g)["__julia_type__"] = string(typeof(f))
    return g
end
=#
