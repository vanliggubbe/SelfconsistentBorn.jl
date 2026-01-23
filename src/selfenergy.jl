struct SelfEnergy{T <: RationalInterpolation, S <: Real}
    dim :: Int
    shift :: S
    Σ :: T
end

evaluate!(y, Σ :: SelfEnergy, x) = lmul!(-1.0im, evaluate!(y, Σ.Σ, x - im * Σ.shift))
(f :: SelfEnergy)(x) = f.Σ(x - im * f.shift) * (-im)

function self_energy_0(oqs :: OQSystem)
    tmp = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 1)
    fill!(tmp, zero(eltype(tmp)))
    # for b in oqs.baths
    #     @inbounds for (j, left) in enumerate(b.cpl_q)
    #         @inbounds for (k, right) in enumerate(b.cpl_c)
    #             mul!(slice(tmp, 1), left' * right, 1.0im, b.R.cnst[j, k], one(eltype(tmp)))
    #         end
    #     end   
    # end
    return SelfEnergy(oqs.dim, 0.0, BarycentricInterpolation([0.0], [1.0], tmp))
end

struct ConstVector{T} <: AbstractVector{T}
    size :: Int
    val :: T
end

Base.size(a :: ConstVector) = (a.size,)
Base.getindex(a :: ConstVector{T}, _) where {T} = a.val

function selfconsistency!(y, oqs :: OQSystem, Σ :: SelfEnergy, ω)
    fill!(y, zero(eltype(y)))
    Id = Diagonal(ConstVector(oqs.dim ^ 2, -one(ComplexF64)))
    ws = LUWs(oqs.ginv)
    
    for b in oqs.baths
        local rR, rK
        #lefts = [Matrix(left') for left in]
        @inbounds for (i, p) in enumerate(b.K.poles)
            if i <= lastdim(b.R.residues)
                rR = slice(b.R.residues, i)
            end
            rK = slice(b.K.residues, i)
            evaluate!(oqs.ginv, Σ, ω - p)
            oqs.ginv .+= oqs.L
            @inbounds for j in axes(oqs.ginv, 1)
                oqs.ginv[j, j] += (p - ω)
            end
            #lmul!(-1.0, ginv)
            F = LU(LAPACK.getrf!(ws, oqs.ginv)...)
            ldiv!(oqs.tmp2, F, Id)
            @inbounds for (j, left) in enumerate(b.cpl_q)
                #mul!(tmp, left', -one(T))
                #rdiv!(tmp, F)
                mul!(oqs.tmp1, left', oqs.tmp2)
                @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                    if i <= lastdim(b.R.residues)
                        mul!(y, oqs.tmp1, right, rR[j, k], one(eltype(y)))
                    end
                    mul!(y, oqs.tmp1, right′, rK[j, k], one(eltype(y)))
                end
            end
        end
        # counterterm
        @inbounds for (k, right) in enumerate(b.cpl_c)
            mul!(oqs.tmp1, right, Id)
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(y, left', oqs.tmp1, -b.R.cnst[j, k], one(eltype(y)))
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
    Ω_cutoff :: Real,
    η :: Real = 0.0;
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)
    @assert (η >= 0.0)
    F! = let oqs = oqs, Σ = Σ; (y, ω) -> lmul!(1.0im, selfconsistency!(y, oqs, Σ, ω + im * η)) end
    F_int = _simple_iteration(oqs, F!, Ω_cutoff; aaa_iter, aaa_eps, aaa_split, nrm)
    return SelfEnergy(oqs.dim, η, F_int)
end