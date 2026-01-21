struct SelfEnergy{T <: RationalInterpolation}
    dim :: Int
    Σ :: T
end

evaluate!(y, Σ :: SelfEnergy, x) = lmul!(-1.0im, evaluate!(y, Σ.Σ, x))
(f :: SelfEnergy)(x) = f.Σ(x) * (-im)

function self_energy_0(oqs :: OQSystem)
    return SelfEnergy(oqs.dim, BarycentricInterpolation([0.0], [1.0], ElasticArray(zeros(oqs.dim ^ 2, oqs.dim ^ 2, 1))))
end

function selfconsistency!(y, oqs :: OQSystem, Σ :: SelfEnergy, ω)
    fill!(y, zero(eltype(y)))
    T = ComplexF64
    ginv = zeros(T, oqs.dim ^ 2, oqs.dim ^ 2)
    tmp = zero(ginv)
    ws = LUWs(ginv)
    L = Matrix(Commutator(oqs.H, -1))
    
    for b in oqs.baths
        local rR, rK
        @inbounds for (i, p) in enumerate(b.K.poles)
            if i <= lastdim(b.R.residues)
                rR = slice(b.R.residues, i)
            end
            rK = slice(b.K.residues, i)
            evaluate!(ginv, Σ, ω - p)
            axpby!(-1.0, L, -1.0, ginv)
            @inbounds for j in axes(ginv, 1)
                ginv[j, j] += ω - p
            end
            F = LU(LAPACK.getrf!(ws, ginv)...)

            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(tmp, left', 1.0)
                rdiv!(tmp, F)
                @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                    if i <= lastdim(b.R.residues)
                        mul!(y, tmp, right, rR[j, k], one(eltype(y)))
                    end
                    mul!(y, tmp, right′, rK[j, k], one(eltype(y)))
                end
            end
        end
        # counterterm
        @inbounds for (j, left) in enumerate(b.cpl_q)
            @inbounds for (k, right) in enumerate(b.cpl_c)
                mul!(y, left' * right, -one(eltype(right)), b.R.cnst[j, k], one(eltype(y)))
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
    _, _, V = svd!(Matrix(Commutator(oqs.H)) + Σ(0.0))
    ret = reshape(V[:, end], oqs.dim, oqs.dim)
    ret ./= tr(ret)
    return ret
end

"""
    simple_iteration(oqs :: OQSystem, Σ :: SelfEnergy, Ω_cutoff :: Real[, which]; kwargs...)

Iterate selfconsistensy equation on self-energy `Σ` for the system `oqs`.
Returns barycentric interpolation on the segment from `-Ω_cutoff` to `Ω_cutoff`.

Keyword arguments:
- `aaa_iter`: number of AAA algorithm iterations, default 20
- `aaa_eps`: accuracy threshold for AAA algorithm, default 1e-9
- `aaa_split`: function which takes iteration index and returns number of splitting points between consecutive support nodes, default `(n -> max(3, 20 - n))`
- `nrm`: norm function, default `LinearAlgebra.norm`
"""
function simple_iteration(
    oqs :: OQSystem,
    Σ :: SelfEnergy,
    Ω_cutoff :: Real;
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)
    F = let oqs = oqs, Σ = Σ; ω -> (selfconsistency(oqs, Σ, ω) * im) end
    F_int = _simple_iteration(oqs, F, Ω_cutoff; aaa_iter, aaa_eps, aaa_split, nrm)
    return SelfEnergy(oqs.dim, F_int)
end