struct BosonicBath{F <: PoleInterpolation, G <: Function, T1 <: Tuple{Vararg{Commutator}}, T2 <: Tuple{Vararg{Commutator}}}
    R :: F
    K :: G

    cpl_c :: T1
    cpl_q :: T2

    function BosonicBath(cpl, R, K)
        #if (lastdim(cpl), lastdim(cpl)) != firstdims(R.residues) || (lastdim(cpl), lastdim(cpl)) != firstdims(K.residues)
        #    throw(DimensionMismatch("Number of coupling operators must coincide with the bath spectral density dimension"))
        #end
        class = Tuple(Commutator(q, 1) for q in cpl)
        quant = Tuple(Commutator(q, -1) for q in cpl)
        new{typeof(R), typeof(K), typeof(class), typeof(quant)}(R, K, class, quant)
    end
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
        n_iter = 100, split = 20, ε = 1e-9,
        weight = (x -> let A = R(x); norm(A - A'); end)
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

"""
    BosonicBath(cpl, J, T, Ω_cutoff)

Creates a BosonicBath object from spectral density `J`, temperature `T`, and list of coupling operators `cpl`.
Function `J` must return a square matrix, it's size must coincide with length of `cpl`.
"""
function BosonicBath(cpl, J :: Function, T :: Real, Ω :: Real, counterterm :: Bool = false)
    R = aaa_real_axis(Ω, x -> J(x) / x)
    K = aaa_real_axis(Ω, x -> J(x) * coth(x / (2 * T)))
    retardize!(R)
    retardize!(K)
    for i in eachindex(R.poles)
        slice(R.residues, i) .*= -im * R.poles[i]
    end
    K.residues .*= -im
    if counterterm
        R.cnst .-= R(0.0)
    end
    #=
    if coth == :aaa
        K = retarded_to_keldysh(R, T, Ω_cutoff)
    elseif coth == :digamma
        K = Keldysh(R, T)
    end
    =#
    return BosonicBath(
        cpl,
        R, K
    )
end

function add_counterterm!(b :: BosonicBath)
    b.R.cnst .= -b.R(0.0)
    nothing
end

#couplings(b :: BosonicBath) = (slice(b.cpl, i) for i in axes(b.cpl, 3))

bose_R(ω, T) = -digamma(ω / (2π * im * T) + 1.0) / (im * π)
bose_A(ω, T) = digamma(-ω / (2π * im * T) + 1.0) / (im * π)

struct Keldysh{S <: Real, TR <: PoleInterpolation} <: Function
    R :: TR
    A :: TR
    R′ :: TR
    A′ :: TR
    T :: S

    function Keldysh(R, T)
        A = PoleInterpolation(
            conj(R.poles), similar(R.residues), Matrix(R.cnst')
        )
        @inbounds for i in axes(A.residues, 3)
            adjoint!(slice(A.residues, i), slice(R.residues, i))
        end

        R′ = PoleInterpolation(R.poles, copy(R.residues), zero(R.cnst))
        for (p, r) in zip(poles(R′), residues(R′))
            lmul!(bose_A(p, T) + (2.0 * T / p), r)
        end
        A′ = PoleInterpolation(A.poles, copy(A.residues), zero(A.cnst))
        for (p, r) in zip(poles(A′), residues(A′))
            lmul!(bose_R(p, T), r)
        end
                
        new{typeof(T), typeof(R)}(R, A, R′, A′, T)
    end
end


function (F :: Keldysh)(ω)
    y = zeros(ComplexF64, firstdims(F.R.residues))
    evaluate!(y, F, ω)
end

function evaluate!(y, F :: Keldysh, ω)
    b = bose_R(ω, F.T)
    evaluate!(y, F.R, ω)
    lmul!(b, y)
    evaluate!(y, F.R′, ω, 1.0, one(eltype(y)))
    evaluate!(y, F.A, ω, -b, one(eltype(y)))
    evaluate!(y, F.A′, ω, 1.0, one(eltype(y)))
    y
end
