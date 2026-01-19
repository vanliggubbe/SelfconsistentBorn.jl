struct PermutationMap{P} <: LinearMap{Bool}
    p :: P

    function PermutationMap(p)
        if !isperm(p)
            throw(ArgumentError("Not a permutation"))
        end
        new{typeof(p)}(p)
    end
end

Base.size(p :: PermutationMap) = (length(p.p), length(p.p))
Base.iterate(p :: PermutationMap) = iterate(p.p)
Base.iterate(p :: PermutationMap, i) = iterate(p.p, i)

MulStyle(::PermutationMap) = FiveArg()

function _unsafe_mul!(y, P :: PermutationMap, x :: AbstractVector)
    y .= @view x[P.p]
    y
end

function _unsafe_mul!(y, P :: AdjointMap{<: Any, <: PermutationMap}, x :: AbstractVector) 
    @views y[P.parent.p] .= x
    y
end

function _unsafe_mul!(y, P :: TransposeMap{<: Any, <: PermutationMap}, x :: AbstractVector) 
    @views y[P.parent.p] .= x
    y
end

function _unsafe_mul!(y, P :: PermutationMap, x :: AbstractVector, α, β) 
    @views axpby!(α, x[P.p], β, y)
    y
end

function _unsafe_mul!(y, P :: AdjointMap{<: Any, <: PermutationMap}, x :: AbstractVector, α, β) 
    @views axpby!(α, x, β, y[P.parent.p])
    y
end

function _unsafe_mul!(y, P :: TransposeMap{<: Any, <: PermutationMap}, x :: AbstractVector, α, β) 
    @views axpby!(α, x, β, y[P.parent.p])
    y
end

"""
    evaluate!(Y, P :: PermutationMap, X :: AbstractMatrix)

Overwrite `Y` with `P * X * P'`
"""
function evaluate!(y, p :: PermutationMap, x :: AbstractMatrix) 
    (@views y[:, p.p] .= x[p.p, :])
    y
end

"""
    evaluate!(Y, P :: PermutationMap, X :: AbstractMatrix, α :: Number, β :: Number)

Overwrite `Y` with  `α * (P * Y * P') + β * Y`
"""
function evaluate!(y, p :: PermutationMap, x :: AbstractMatrix, α :: Number, β :: Number)
    (@views axpby!(α, x[p.p, :], β, y[:, p.p]))
    y
end

"""
    evaluate!(Y, P :: ComposedFunction{typeof(conj), <: PermutationMap}, X :: AbstractMatrix)

Overwrite `Y` with `conj(P * X * P')`
"""
function evaluate!(y, p :: ComposedFunction{typeof(conj), <: PermutationMap}, x :: AbstractMatrix)
    evaluate!(y, p.inner, x)
    conj!(y)
end

"""
    evaluate!(Y, P :: ComposedFunction{typeof(conj), <: PermutationMap}, X :: AbstractMatrix, α :: Number, β :: Number)

Overwrite `Y` with `α * conj(P * X * P') + β * Y`
"""
function evaluate!(y, p :: ComposedFunction{typeof(conj), <: PermutationMap}, x :: AbstractMatrix, α :: Number, β :: Number)
    conj!(y)
    evaluate!(y, p.inner, x, conj(α), conj(β))
    conj!(y)
end
