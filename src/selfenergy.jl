struct SelfEnergy{T <: RationalInterpolation{3}, S <: Real, LIOU <: AbstractMatrix}
    dim :: Int
    shift :: S
    Σ :: T
    L :: LIOU
    idxs :: Vector{Int}
    blocks :: Vector{UnitRange{Int}}
end


function block_structure!(blocks :: IntDisjointSet, f :: PoleInterpolation{3}; ε :: Real = 1e-12)
    block_structure!(blocks, f.cnst; ε)
    block_structure!(blocks, f.residues; ε)
end

function block_structure!(blocks :: IntDisjointSet, f :: BarycentricInterpolation{3}; ε :: Real = 1e-12)
    @inbounds for n in axes(f.values, 3)
        block_structure!(blocks, slice(f.values, n); ε)
    end
end

function block_structure!(blocks :: IntDisjointSet, f :: SymmetricBarycentricInterpolation{3}; ε :: Real = 1e-12)
    A = Matrix{eltype(f.values)}(undef, firstdims(f.values))
    @inbounds for n in axes(f.values, 3)
        block_structure!(blocks, slice(f.values, n); ε)
        evaluate!(A, f.sym_v, slice(f.values, n))
        block_structure!(blocks, A)
    end
end

function Base.permute!(f :: Union{<: BarycentricInterpolation, SymmetricBarycentricInterpolation}, perm :: Vector{Int})
    tmp = Matrix{eltype(f.values)}(undef, firstdims(f.values))
    for i in axes(f.values, 3)
        tmp .= slice(f.values, i)[perm, perm]
        slice(f.values, i) .= tmp
    end
end

function Base.permute!(f :: PoleInterpolation, perm :: Vector{Int})
    tmp = Matrix{eltype(f.residues)}(undef, firstdims(f.residues))
    tmp .= f.cnst[perm, perm]
    f.cnst .= tmp
    for i in axes(f.residues, ndims(f.residues))
        tmp .= slice(f.residues, i)[perm, perm]
        slice(f.residues, i) .= tmp
    end
end


function SelfEnergy(dim :: Int, η :: Real, Σ :: RationalInterpolation{3}, L :: AbstractMatrix; ε :: Real = 1e-12)
    # identify matrix blocks
    # we do it each time new self-energy is constructed
    blocks = IntDisjointSet(dim ^ 2)
    block_structure!(blocks, L; ε)
    block_structure!(blocks, Σ; ε)
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
    permute!(Σ, idxs)
    if Σ isa SymmetricBarycentricInterpolation
        rmul!(lmul!(PermutationMap(idxs), Σ.sym_v.inner), PermutationMap(invperm(idxs)))
    end
    SelfEnergy(dim, η, Σ, L[idxs, idxs], idxs, qqq)
end

evaluate!(y, Σ :: SelfEnergy, x) = evaluate!(y, Σ.Σ, x - im * Σ.shift)
evaluate!(y, Σ :: SelfEnergy, x, :: Val{true}) = evaluate!(view(y, Σ.idxs, Σ.idxs), Σ.Σ, x - im * Σ.shift)

function (f :: SelfEnergy)(x)
    y = Matrix{ComplexF64}(undef, f.dim ^ 2, f.dim ^ 2)
    evaluate!(y, f, x, Val{true}())
    y
end

function self_energy_0(oqs :: OQSystem)
    tmp = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 1)
    fill!(tmp, zero(eltype(tmp)))
    return SelfEnergy(oqs.dim, 0.0, BarycentricInterpolation([0.0], [1.0], tmp), oqs.L)
end

function min_q(β :: AbstractVector{<: Real}, ζ :: AbstractVector{<: Real}, σ :: Real; n_ternary :: Integer = 30)
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
        _, neg = min_q(β, ζ, σ_r)
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
        _, neg = min_q(β, ζ, σ)
        if neg
            σ_r = σ
        else
            σ_l = σ
        end
    end
    return σ_r
end

function counterterm!(y, oqs :: OQSystem)
    # counterterm
    for b in oqs.baths
        @inbounds for (k, right) in enumerate(b.cpl_c)
            mul!(oqs.tmp1, right, I(size(right, 1)))
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(y, left', oqs.tmp1, b.R.cnst[j, k], one(eltype(y)))
            end
        end    
    end 
    return y
end

function selfconsistency!(y, oqs :: OQSystem, Σ :: SelfEnergy, ω)
    fill!(y, zero(eltype(y)))
    counterterm!(y, oqs)
    ws = LUWs(oqs.ginv)
    Is = [I(length(x)) for x in Σ.blocks]
    
    for b in oqs.baths
        # calculate separately retarded and Keldysh poles
        for (S, cpl_l, cpl_r) in [(b.R, b.cpl_c, b.cpl_q), (b.K, b.cpl_q, b.cpl_q)]
            @inbounds for (i, p) in enumerate(S.poles)
                rS = slice(S.residues, i)
                evaluate!(oqs.ginv, Σ, ω - p)
                oqs.ginv .+= Σ.L
                @inbounds for j in axes(oqs.ginv, 1)
                    oqs.ginv[j, j] += (p - ω)
                end
                fill!(oqs.tmp1, zero(eltype(oqs.tmp1)))
                for (II, rg) in zip(Is, Σ.blocks)
                    if length(rg) > 1
                        F = LU(LAPACK.getrf!(ws, view(oqs.ginv, rg, rg))...)
                        ldiv!(view(oqs.tmp1, rg, rg), F, II)
                    else
                        oqs.tmp1[first(rg), first(rg)] = inv(oqs.ginv[first(rg), first(rg)])
                    end
                end
                @views oqs.tmp2[Σ.idxs, Σ.idxs] .= oqs.tmp1

                @inbounds for (j, left) in enumerate(cpl_l)
                    mul!(oqs.tmp1, left', oqs.tmp2)
                    @inbounds for (k, right) in enumerate(cpl_r)
                        # minus sign because ginv has a wrong sign
                        mul!(y, oqs.tmp1, right, -rS[j, k], one(eltype(y)))
                    end
                end
            end
        end
   end
   return y
end

function selfconsistency(oqs :: OQSystem, Σ :: SelfEnergy, ω)
    y = Matrix{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2)
    selfconsistency!(y, oqs, Σ, ω)
end

"""
    steady_state(oqs :: OQSystem, Σ :: SelfEnergy)

Returns steady state density matrix given the open quantum system and converged self-energy
"""
function steady_state(oqs :: OQSystem, Σ :: SelfEnergy)
    _, _, V = svd(oqs.L + Σ(0.0))
    ret = reshape(V[:, end], oqs.dim, oqs.dim)
    ret ./= tr(ret)
    return ret
end

"""
    simple_iteration(oqs :: OQSystem, Σ :: SelfEnergy, Ω_cutoff :: Real[, η :: Real = 0.0]; kwargs...)

Iterate selfconsistensy equation on self-energy `Σ` for the system `oqs`.
Returns barycentric interpolation on the segment from `-Ω_cutoff + im * η` to `Ω_cutoff + im * η`.
Imaginary displacement `η` must be non-negative.

Keyword arguments:
- `aaa_iter`: number of AAA algorithm iterations, default 20
- `aaa_eps`: accuracy threshold for AAA algorithm, default 1e-9
- `aaa_split`: function which takes iteration index and returns number of splitting points between consecutive support nodes, default `(n -> max(3, 20 - n))`
- `nrm`: norm function, default `LinearAlgebra.norm`
"""
function simple_iteration(
    oqs :: OQSystem,
    Σ :: SelfEnergy,
    Ω :: Real,
    η :: Real = simple_iteration_domain(oqs);
    simple_iter :: Int = 5,
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
)
    @argcheck (η >= 0.0)

    # construct mutating functions and symmetry operator
    pm = PermutationMap(copy(oqs.perm.p))
    Σ′ = Σ
    sym!(y, x) = lmul!(-1, conj!(evaluate!(y, pm, x)))
    for _ in 1 : simple_iter
        F!(y, ω) = selfconsistency!(y, oqs, Σ′, ω + im * η)

        F_int = aaa_real_axis(
            PoleInterpolation,
            Ω, F!, sym!;
            aaa_iter, aaa_eps = (aaa_eps * oqs.dim ^ 2), aaa_split,
            fun_mut = true, fun_shape = (oqs.dim ^ 2, oqs.dim ^ 2), fun_type = ComplexF64,
            sym_mut = true
        )
        counterterm!(F_int.cnst, oqs)
        Σ′ = SelfEnergy(oqs.dim, η, F_int, oqs.L)
    end
    return Σ′
end

# struct SelfEnergy{T <: RationalInterpolation{3}, S <: Real, LIOU <: AbstractMatrix}
#     dim :: Int
#     shift :: S
#     Σ :: T
#     L :: LIOU
#     idxs :: Vector{Int}
#     blocks :: Vector{UnitRange{Int}}
# end

function Base.write(parent :: Union{File, Group}, name :: AbstractString, f :: SelfEnergy)
    g = create_group(parent, name)
    g["dim"] = f.dim
    g["shift"] = f.shift
    g["self-energy"] = f.Σ
    g["liou"] = f.L
    g["idxs"] = f.idxs

    attributes(g)["__julia_type__"] = string(typeof(f))
    return g
end

function Base.read(parent :: Union{File, Group}, name :: AbstractString, :: Type{<: SelfEnergy})
    g = parent[name]
    idxs = read(g, "idxs")
    cur = 1
    blocks = UnitRange{Int}[]
    for i in eachindex(idxs[begin : end - 1])
        if i == lastindex(idxs) || idxs[i + 1] < idxs[i]
            push!(blocks, cur : i)
            cur = i + 1
        end
    end
    SelfEnergy(
        read(g, "dim"),
        read(g, "shift"),
        read(g, "self-energy", SymmetricBarycentricInterpolation),
        read(g, "liou"),
        idxs, blocks
    )
end
