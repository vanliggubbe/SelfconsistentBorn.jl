struct Commutator{T, M <: Union{AbstractMatrix{T}, LinearMap{T}}} <: LinearMap{T}
    A :: M
    n :: Int
    sign :: T

    function Commutator(A, s)
        if size(A, 1) != size(A, 2)
            throw(DimensionMismatch("Matrix must be square"))
        end
        new{eltype(A), typeof(A)}(A, size(A, 1), convert(eltype(A), sign(s)))
    end
end

Commutator(A) = Commutator(A, -1)

ishermitian(A :: Commutator) = ishermitian(A.A)
issymmetric(A :: Commutator) = issymmetric(A.A)
#Base.adjoint(A :: Commutator) = Commutator(A.A', A.sign)
#Base.transpose(A :: Commutator) = Commutator(transpose(A.A), A.sign)
Base.size(A :: Commutator) = (A.n ^ 2, A.n ^ 2)

MulStyle(::Commutator) = FiveArg()

function _unsafe_mul!(y, A :: Commutator, x :: AbstractVector)
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, A.A, X, one(eltype(X)), zero(eltype(y)))
    mul!(Y, X, A.A, A.sign, one(eltype(y)))
    return y
end

function _unsafe_mul!(y, B :: AdjointMap{<: Any, <: Commutator}, x :: AbstractVector)
    A = B.lmap
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, A.A', X, one(eltype(X)), zero(eltype(y)))
    mul!(Y, X, A.A', A.sign, one(eltype(y)))
    return y
end

function _unsafe_mul!(y, B :: TransposeMap{<: Any, <: Commutator}, x :: AbstractVector)
    A = B.lmap
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, transpose(A.A), X, one(eltype(X)), zero(eltype(y)))
    mul!(Y, X, transpose(A.A), A.sign, one(eltype(y)))
    return y
end

function _unsafe_mul!(y, A :: Commutator, x :: AbstractVector, α, β)
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, A.A, X, α, β)
    mul!(Y, X, A.A, A.sign * α, one(eltype(Y)))
    return y
end

function _unsafe_mul!(y, B :: AdjointMap{<: Any, <: Commutator}, x :: AbstractVector, α, β)
    A = B.lmap
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, A.A', X, α, β)
    mul!(Y, X, A.A', A.sign * α, one(eltype(Y)))
    return y
end

function _unsafe_mul!(y, B :: TransposeMap{<: Any, <: Commutator}, x :: AbstractVector, α, β)
    A = B.lmap
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, transpose(A.A), X, α, β)
    mul!(Y, X, transpose(A.A), A.sign * α, one(eltype(Y)))
    return y
end