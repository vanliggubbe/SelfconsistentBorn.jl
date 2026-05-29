@inline @generated function divtype(:: Type{T}, :: Type{S}) where {T <: Number, S <: Number}
    ret = typeof(one(T) / one(S))
    return :($ret)
end

@inline @generated function divtype(:: T, :: S) where {T <: Number, S <: Number}
    ret = typeof(one(T) / one(S))
    return :($ret)
end

@inline slice(a :: AbstractArray{T, 1}, i :: Integer) where {T} = a[i]
@inline slice(a :: AbstractArray{T, N}, i) where {T, N} = selectdim(a, N, i)
@inline lastdim(a :: AbstractArray{T, N}) where {T, N} = size(a, N)
@inline firstdims(a :: AbstractArray{T, N}) where {T, N} = size(a)[1 : N - 1]

@inline flatten(x :: Number) = x
@inline flatten(x :: Array) = vec(x)

block_structure!(:: IntDisjointSet, :: AbstractVector; ε :: Real = 1e-12) = error("Not applicable")

function block_structure!(blocks :: IntDisjointSet, a :: AbstractMatrix; ε :: Real = 1e-12)
    @inbounds for i in axes(a, 2)
        @inbounds for j in axes(a, 1)
            if abs(a[j, i]) > ε
                union!(blocks, i, j)
            end
        end
    end
end

function block_structure!(block :: IntDisjointSet, a :: AbstractArray; ε :: Real = 1e-12)
    @inbounds for i in axes(a, ndims(a))
        block_structure!(block, slice(a, i); ε)
    end
end

intermediate(Ω :: Real, a :: Real, b :: Real, n :: Int) = Ω * tan.(
    LinRange(
        angle((im * Ω - a) / (im * Ω + a)),
        isfinite(b) ? angle((im * Ω - b) / (im * Ω + b)) : π,
        n + 2
    )[2 : end - 1] / 2.0
)

intermediate(a :: Complex, b :: Complex, n :: Int) = exp.(
    LinRange(log(a), log(b), n + 2)[2 : end - 1]
)

#evaluate!(y, f :: Function, x :: Number) = (y .= f(x))
#evaluate!(y, f :: Function, x :: AbstractArray) = map!(f, y, x)

minusconj(x) = -conj(x)

hermitian_part(a) = (a + a') / 2
antihermitian_part(a) = (a - a') / 2
