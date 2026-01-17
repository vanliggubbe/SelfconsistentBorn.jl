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

Base.adjoint(A :: Commutator) = Commutator(A.A', A.sign)
Base.transpose(A :: Commutator) = Commutator(transpose(A.A), A.sign)
Base.size(A :: Commutator) = (A.n ^ 2, A.n ^ 2)

function _unsafe_mul!(y, A :: Commutator, x :: AbstractVector)
    Y = reshape(y, A.n, A.n)
    X = reshape(x, A.n, A.n)
    mul!(Y, A.A, X)
    mul!(Y, X, A.A, A.sign, one(eltype(y)))
    return y
end