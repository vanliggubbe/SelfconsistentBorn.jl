struct SelfEnergySymmetry{P <: PermutationMap} <: Function
    pm :: P
end

evaluate!(y, s :: SelfEnergySymmetry, x) = lmul!(-one(eltype(y)), conj!(evaluate!(y, s.pm, x)))
# y := s(x) * a + b * y
evaluate!(y, s :: SelfEnergySymmetry, x, a, b) = lmul!(-one(eltype(y)), conj!(evaluate!(conj!(y), s.pm, x, conj(a), -conj(b))))
(f :: SelfEnergySymmetry)(x) = (evaluate!(similar(x), f, x))
(f :: SelfEnergySymmetry)(y, x) = evaluate!(y, f, x)

struct OQSystem{
    HAMILTONIAN <: Union{AbstractMatrix, LinearMap},
    BATHS,
    SELFENERGY <: SymmetricBarycentricInterpolation{3}
}
    dim         :: Int
    hamiltonian :: HAMILTONIAN
    baths       :: BATHS
    selfenergy  :: SELFENERGY

    idxs        :: Vector{Int}
    blocks      :: Vector{UnitRange{Int}}
    liou        :: Matrix{ComplexF64}

    function OQSystem(H, baths; ε :: Real = 1e-12)
        @argcheck ε > zero(eltype(ε))
        if size(H, 1) != size(H, 2)
            throw(DimensionMismatch("Hamiltonian must be square"))
        end
        dim = size(H, 1)
        for b in baths
            if !(b isa BosonicBath)
                throw(ArgumentError("Each element of `bath` must have a type of bosonic bath"))
            end
            for (c, q) in zip(b.cpl_c, b.cpl_q)
                if size(c) != (dim ^ 2, dim ^ 2) || size(q) != (dim ^ 2, dim ^ 2)
                    throw(DimensionMismatch("Dimensions of the Hamiltonian and the coupling operators must coincide"))
                end
            end
        end
        # liouvillian
        L = Matrix{ComplexF64}(undef, dim ^ 2, dim ^ 2)
        mul!(L, Commutator(H), one(eltype(L)))

        perm = collect(vec(transpose(reshape(1 : dim ^ 2, dim, dim))))

        # block structure of the OQS liouvillian
        blocks = IntDisjointSet(dim ^ 2)
        block_structure!(blocks, L; ε)
        for b in baths
            for cpl_q in b.cpl_q
                for (cpl_c′, cpl_q′) in zip(b.cpl_c, b.cpl_q)
                    block_structure!(blocks, Matrix(cpl_q * cpl_c′); ε)
                    block_structure!(blocks, Matrix(cpl_q * cpl_q′); ε)
                end
            end
        end

        bl = Vector{Int}[]
        blmap = zeros(Int, dim ^ 2)
        for i in 1 : dim ^ 2
            r = find_root!(blocks, i)
            if blmap[r] == 0
                push!(bl, Int[])
                blmap[r] = length(bl)
            end
            push!(bl[blmap[r]], i)
        end
        cur = 0
        qqq = UnitRange{Int}[]
        for block in bl
            push!(qqq, cur + 1 : cur + length(block))
            cur += length(block)
        end
        idxs = reduce(vcat, bl)

        # construct zero-valued self-energy
        sym_Σ = SelfEnergySymmetry(PermutationMap(perm))
        Σ = SymmetricBarycentricInterpolation(
            ComplexF64[],
            ones(ComplexF64, 1),
            zeros(ComplexF64, dim ^ 2, dim ^ 2, 1),
            minusconj,
            minusconj,
            sym_Σ
        )

        new{typeof(H), typeof(baths), typeof(Σ)}(
            dim, H, baths, Σ, invperm(idxs), qqq, L
        )
    end
end

@inline dim(oqs :: OQSystem) = oqs.dim

self_energy!(y, oqs :: OQSystem, x) = evaluate!(y, oqs.selfenergy, x)

function self_energy(oqs :: OQSystem, x)
    y = Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2)
    self_energy!(y, oqs, x)
    y
end

function _min_q(β :: AbstractVector{<: Real}, ζ :: AbstractVector{<: Real}, σ :: Real; n_ternary :: Integer = 30)
    ζ′ = ζ .+ σ
    X = sum(β ./ ζ′)
    q_r = 1.0 - X / minimum(ζ′)
    if q_r <= 0.0
        return 0.0, false
    end

    # ternary search
    q_l = 0.0
    for _ in 1 : n_ternary
        q_l′ = q_l + (q_r - q_l) / 3.0
        q_r′ = q_r - (q_r - q_l) / 3.0
        f_l′ = sum(β ./ (ζ′ .- X / (1 - q_l′)) .^ 2) - q_l′
        f_r′ = sum(β ./ (ζ′ .- X / (1 - q_r′)) .^ 2) - q_r′
        if f_l′ < f_r′
            q_r = q_r′
        else
            q_l = q_l′
        end
    end
    f_l = sum(β ./ (ζ′ .- X / (1 - q_l)) .^ 2) - q_l
    f_r = sum(β ./ (ζ′ .- X / (1 - q_r)) .^ 2) - q_r
    if f_l < f_r
        return q_l, (f_l < 0.0)
    else
        return q_r, (f_r < 0.0)
    end
end

# calculate domain that simple iteration converges
function simple_iteration_domain(oqs :: OQSystem)
    β = Float64[]
    ζ = Float64[]
    # get bounds for coupling operators
    for b in oqs.baths
        norm_c = norm_bound.(b.cpl_c)
        norm_q = norm_bound.(b.cpl_q)
        @inbounds for i in eachindex(b.R.poles)
            push!(ζ, -imag(b.R.poles[i]))
            rR = slice(b.R.residues, i)
            tmp = 0.0
            for (j, a) in enumerate(norm_c)
                for (k, b) in enumerate(norm_q)
                    tmp += a * b * abs(rR[j, k])
                end
            end
            push!(β, tmp)
        end
        @inbounds for i in eachindex(b.K.poles)
            push!(ζ, -imag(b.K.poles[i]))
            rK = slice(b.K.residues, i)
            tmp = 0.0
            for (j, a) in enumerate(norm_q)
                for (k, b) in enumerate(norm_q)
                    tmp += a * b * abs(rK[j, k])
                end
            end
            push!(β, tmp)
        end
    end

    σ_r = 0.0
    while true
        _, neg = _min_q(β, ζ, σ_r)
        if neg
            break
        end
        σ_r = 2.0 * σ_r + 1.0
    end
    if iszero(σ_r)
        return 0.0
    end
    σ_l = (σ_r - 1.0) / 2.0
    # binary search
    for _ in 1 : 20
        σ = (σ_l + σ_r) / 2.0
        _, neg = _min_q(β, ζ, σ)
        if neg
            σ_r = σ
        else
            σ_l = σ
        end
    end
    return σ_r
end

function counterterm!(y, oqs :: OQSystem; ws = Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2))
    # counterterm
    for b in oqs.baths
        @inbounds for (k, right) in enumerate(b.cpl_c)
            mul!(ws, right, I(size(right, 1)))
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(y, left', ws, b.R.cnst[j, k], one(eltype(y)))
            end
        end    
    end 
    return y
end

function selfconsistency!(
    y, oqs :: OQSystem, ω;
    ws = (
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2),
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2),
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2)
    )
)
    fill!(y, zero(eltype(y)))
    counterterm!(y, oqs; ws = ws[2])
    luw_ws = LUWs(ws[1])
    Is = [I(length(x)) for x in oqs.blocks]
    ginv, tmp1, tmp2 = ws
    
    for b in oqs.baths
        # calculate separately retarded and Keldysh poles
        for (S, cpl_l, cpl_r) in [(b.R, b.cpl_q, b.cpl_c), (b.K, b.cpl_q, b.cpl_q)]
            @inbounds for (i, p) in enumerate(S.poles)
                rS = slice(S.residues, i)
                @views self_energy!(ginv[oqs.idxs, oqs.idxs], oqs, ω - p)
                @views ginv[oqs.idxs, oqs.idxs] .+= oqs.liou
                @inbounds for j in axes(ginv, 1)
                    ginv[j, j] += (p - ω)
                end
                fill!(tmp1, zero(eltype(tmp1)))
                for (II, rg) in zip(Is, oqs.blocks)
                    if length(rg) > 1
                        F = LU(LAPACK.getrf!(luw_ws, view(ginv, rg, rg))...)
                        ldiv!(view(tmp1, rg, rg), F, II)
                    else
                        tmp1[first(rg), first(rg)] = inv(ginv[first(rg), first(rg)])
                    end
                end
                @views tmp2 .= tmp1[oqs.idxs, oqs.idxs]

                @inbounds for (j, left) in enumerate(cpl_l)
                    mul!(tmp1, left', tmp2)
                    @inbounds for (k, right) in enumerate(cpl_r)
                        # minus sign because ginv has a wrong sign
                        mul!(y, tmp1, right, -rS[j, k], one(eltype(y)))
                    end
                end
            end
        end
   end
   return y
end

function selfconsistency(oqs :: OQSystem, ω)
    y = Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2)
    selfconsistency!(y, oqs, ω)
end

selfconsistency_check(oqs :: OQSystem, ω; nrm = norm) = nrm(self_energy(oqs, ω) - selfconsistency(oqs, ω))

"""
    steady_state(oqs :: OQSystem)

Returns steady state density matrix given the open quantum system
"""
function steady_state(oqs :: OQSystem)
    _, _, V = svd(oqs.liou + self_energy(oqs, 0.0))
    ret = reshape(V[:, end], dim(oqs), dim(oqs))
    ret ./= tr(ret)
    return ret
end

"""
    simple_iteration(oqs :: OQSystem, Ω :: Real[, η :: Real = simple_iteration_domain(oqs)]; kwargs...)

Iterate selfconsistensy for the open quantum system `oqs`.
Imaginary displacement `η` must be non-negative.

Keyword arguments:
- `simple_iter`: number of simple iterations
- `aaa_iter`: number of AAA algorithm iterations, default 20
- `aaa_eps`: accuracy threshold for AAA algorithm, default 1e-9
- `aaa_split`: function which takes iteration index and returns number of splitting points between consecutive support nodes, default `(n -> max(3, 20 - n))`
"""
function simple_iteration!(
    oqs :: OQSystem,
    Ω   :: Real,
    η   :: Real = simple_iteration_domain(oqs);
    simple_iter :: Int = 5,
    aaa_iter    :: Int = 20,
    aaa_eps     :: Real = 1e-9,
    aaa_split   :: Function = (n -> max(3, 20 - n)),
)
    @argcheck (η >= 0.0)

    # workspace for self-energy evaluation
    ws = (
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2),
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2),
        Matrix{ComplexF64}(undef, dim(oqs) ^ 2, dim(oqs) ^ 2)
    )
    #sym = SelfEnergySymmetry(PermutationMap(vec(transpose(reshape(1 : dim(oqs) ^ 2, dim(oqs), dim(oqs))))))
    for _ in 1 : simple_iter
        F!(y, ω) = selfconsistency!(y, oqs, ω + im * η; ws)

        nodes, weights, values = aaa_real_axis(
            SymmetricBarycentricInterpolation,
            Ω, F!, oqs.selfenergy.sym_v;
            aaa_iter, aaa_eps = (aaa_eps * dim(oqs) ^ 2), aaa_split,
            fun_mut = true, fun_shape = (dim(oqs) ^ 2, dim(oqs) ^ 2), fun_type = ComplexF64,
            sym_mut = true
        )
        resize!(oqs.selfenergy.nodes,   length(nodes))
        resize!(oqs.selfenergy.weights, length(weights))
        resize!(oqs.selfenergy.values,  size(values))
        resize!(oqs.selfenergy.perm,    length(nodes))
        oqs.selfenergy.nodes    .= nodes .+ im * η
        oqs.selfenergy.weights  .= weights
        oqs.selfenergy.values   .= values
        oqs.selfenergy.perm     .= collect(1 : length(nodes))
        sortperm!(oqs.selfenergy.perm, abs.(weights[1 : end - 1]))
    end
    oqs
end


