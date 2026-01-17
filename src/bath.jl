struct BosonicBath{T <: Number, F <: PoleInterpolation, G <: PoleInterpolation}
    cpl :: ElasticArray{T, 3, 2, Vector{T}}
    R :: F
    K :: G
end

function retarded_to_keldysh(R :: PoleInterpolation{3}, T :: Real, Ω_cutoff :: Real; ε :: Real = 1e-15)
    # construct advanced R - A
    residues = ElasticArray{promote_type(eltype(R.residues), ComplexF64)}(undef, size(R.residues, 1), size(R.residues, 2), size(R.residues, 3) * 2)
    slice(residues, 1 : size(R.residues, 3)) .= R.residues
    for i in axes(R.residues, 3)
        adjoint!(slice(residues, i + lastdim(R.residues)), slice(R.residues, i))
        lmul!(-one(eltype(residues)), slice(residues, i + lastdim(R.residues)))
    end
    K = PoleInterpolation(
        [R.poles; conj(R.poles)],
        residues,
        R.cnst - R.cnst'
    )

    # TODO pass parameters for AAA call
    bose = aaa_mat_odd(
        0.125 * T, Ω_cutoff,
        let T = T; x -> (coth(x / (2.0 * T)) - 2.0 * T / x) end;
        n_iter = 100, split = 20, ε = 1e-9
    ) |> PoleInterpolation
    cotanh(x) = bose(x) + 2.0 * T / x

    # multiply (R - A) by cotangent
    residues′ = Array{promote_type(eltype(residues), eltype(bose.poles), eltype(bose.residues), ComplexF64)}(
        undef, size(residues, 1), size(residues, 2), length(bose.poles)
    )
    @inbounds for i in eachindex(bose.poles)
        evaluate!(slice(residues′, i), K, bose.poles[i])
        lmul!(bose.residues[i], slice(residues′, i))
    end
    @inbounds for i in eachindex(K.poles)
        lmul!(cotanh(K.poles[i]), slice(K.residues, i))
    end
    K.cnst .*= bose.cnst
    append!(residues, residues′)
    append!(K.poles, bose.poles)

    # crop poles with advanced causality
    return retardize!(K)
end

function BosonicBath(cpl, J :: Function, T :: Real, Ω_cutoff :: Real)
    sz = size(first(cpl))
    S = eltype(first(cpl))
    R = let tmp = aaa_mat_odd(0.125 * T, Ω_cutoff, J)
        PoleInterpolation(
            OddBarycentricInterpolation(
                tmp.nodes,
                tmp.weights,
                -im * tmp.values
            )
        )
    end
    retardize!(R)
    K = retarded_to_keldysh(R, T, Ω_cutoff)
    return BosonicBath(
        reduce(
            append!, cpl;
            init = ElasticArray{S}(undef, sz..., 0)
        ),
        R, K
    )
end

couplings(b :: BosonicBath) = (slice(b.cpl, i) for i in axes(b.cpl, 3))