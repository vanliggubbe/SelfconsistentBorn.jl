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
