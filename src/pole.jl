mutable struct PoleInterpolation{N, T, A <: ElasticArray{T, N}, C, P <: Number} <: RationalInterpolation{N}
    poles :: Vector{P}
    residues :: A
    cnst :: C
end

function PoleInterpolation(poles, residues :: AbstractVector{Array{T, N}}, cnst :: Array{S, N})  where {T, N, S}
    if length(poles) != length(residues)
        throw(DimensionMismatch("Number of poles must coincide with number of residues"))
    end
    sz = size(cnst)
    ress = reduce(
        append!, residues;
        init = ElasticArray{T}(undef, (sz..., 0))
    )
    return PoleInterpolation(poles isa Vector ? poles : collect(poles), ress, cnst)
end

function PoleInterpolation(poles, residues :: AbstractVector, cnst)
    if length(poles) != length(residues)
        throw(DimensionMismatch("Number of poles must coincide with number of residues"))
    end
    return PoleInterpolation(
        poles isa Vector ? poles : collect(poles),
        residues isa ElasticVector ? residues : ElasticVector(residues),
        cnst
    )
end

PoleInterpolation(poles, residues :: AbstractVector{<: Array}) = PoleInterpolation(poles, residues, zero(first(residues)))
PoleInterpolation(poles, residues :: AbstractVector{T}) where {T} = PoleInterpolation(poles, residues, zero(T))

evaluate!(_, F :: PoleInterpolation{1}, __) = error("Not applicable")

function evaluate!(y, F :: PoleInterpolation, x :: Number)
    y .= F.cnst
    @inbounds for i in eachindex(F.poles)
        axpy!(inv(x - F.poles[i]), slice(F.residues, i), y)
    end
    return y
end

function evaluate!(y, F :: PoleInterpolation, x :: Number, α :: Number, β :: Number)
    axpby!(α, F.cnst, β, y)
    @inbounds for i in eachindex(F.poles)
        axpy!(α / (x - F.poles[i]), slice(F.residues, i), y)
    end
    return y
end

function evaluate_diff!(y, F :: PoleInterpolation, x)
    fill!(y, eltype(y))
    @inbounds for i in eachindex(F.poles)
        axpy!(-inv(x - F.poles[i]) ^ 2, slice(F.residues, i), y)
    end
    return y
end

function (F :: PoleInterpolation{1})(x)
    ret = convert(
        promote_type(
            divtype(eltype(F.residues), promote_type(typeof(x), eltype(F.poles))),
            typeof(F.cnst)
        ),
        F.cnst
    )
    
    @inbounds @simd for i in eachindex(F.poles, F.residues)
        ret += F.residues[i] / (x - F.poles[i])
    end
    return ret
end

function (F :: PoleInterpolation)(x)
    # Float64 because divide
    T = promote_type(
        divtype(
            promote_type(eltype(F.residues)),
            promote_type(typeof(x), eltype(F.poles))
        ),
        eltype(F.cnst)
    )
    y = Array{T}(undef, size(F.cnst))
    y .= F.cnst
    return evaluate!(y, F, x)
end

function retardize!(F :: PoleInterpolation; ε :: Real = 1e-15)
    mask = imag(F.poles) .< -ε
    len = sum(mask)
    F.poles[1 : len] .= F.poles[mask]
    slice(F.residues, 1 : len) .= slice(F.residues, mask)
    resize!(F.poles, len)
    resize!(F.residues, (size(F.cnst)..., len))
    if F.cnst isa Number
        F.cnst /= 2.0
    else
        F.cnst ./= 2.0
    end
    return F
end

# TODO rewrite more efficiently
function PoleInterpolation(F :: OddBarycentricInterpolation)
    # get poles
    B = Matrix(1.0I, length(F.weights) + 1, length(F.weights) + 1)
    B[1, 1] = 0.0
    E = zeros(promote_type(eltype(F.weights), eltype(F.nodes)), size(B))
    E[1, 2 : end] .= sqrt.(abs.(F.weights))
    E[2 : end, 1] .= F.weights ./ E[1, 2 : end]
    E[2 : end, 2 : end] = Diagonal(F.nodes .^ 2)
    pol2 = filter(isfinite, eigvals(E, B))
    pol = sqrt.(Complex.(pol2))
    pol = [pol; -pol]
    pol2 = [pol2; pol2]

    # get residues
    C_pol = inv.(
        (pol2) * ones(length(F.nodes))' -
        ones(length(pol)) * transpose(F.nodes .^ 2)
    )
    N_pol = C_pol * [2.0 * slice(F.values, i) * F.weights[i] * F.nodes[i] for i in 1 : length(F.perm)]
    Ddiff_pol = 2.0 * C_pol * F.weights - 4.0 * (pol2) .* ((C_pol .^ 2) * F.weights)
    N_pol ./= Ddiff_pol

    return PoleInterpolation(
        pol,
        N_pol,
        zero(first(N_pol))
    )
end

function PoleInterpolation(F :: SymmetricBarycentricInterpolation)
    # get poles
    pol = poles(F)

    # get residues
    C_pol = inv.(
        pol * ones(length(F.nodes) * 2)' -
        ones(length(pol)) * transpose([F.nodes; -F.nodes])
    )
    ws = [
        F.weights;
        F.sym_w.(F.weights)
    ]
    vs = cat(F.values, F.values; dims = ndims(F.values))
    for i in axes(F.values, ndims(F.values))
        evaluate!(slice(vs, i + size(F.values, 3)), F.sym_v, slice(F.values, i))
    end
    N_pol = C_pol * [slice(vs, i) * ws[i] for  i in eachindex(ws)]
    Ddiff_pol = (-C_pol .^ 2) * ws
    N_pol ./= Ddiff_pol

    return PoleInterpolation(
        pol,
        N_pol,
        zero(first(N_pol))
    )
end

poles(f :: PoleInterpolation) = f.poles
residues(f :: PoleInterpolation{1}) = f.residues
residues(f :: PoleInterpolation) = (slice(f.residues, i) for i in axes(f.residues, ndims(f.residues)))

function Base.write(parent :: Union{File, Group}, name :: AbstractString, f :: PoleInterpolation)
    g = create_group(parent, name)
    g["poles"] = f.poles
    g["residues"] = f.residues
    g["cnst"] = f.cnst

    attributes(g)["__julia_type__"] = string(typeof(f))
    return g
end

function Base.read(parent :: Union{File, Group}, name :: AbstractString, :: Type{<: PoleInterpolation})
    g = parent[name]
    return PoleInterpolation(
        read(g, "poles"),
        ElasticArray(read(g, "residues")),
        read(g, "cnst")
    )
end
