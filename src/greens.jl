struct GreensFunction{F <: RationalInterpolation{3}}
    dim :: Int
    g :: F
end

function evaluate!(y, F :: GreensFunction, x)
    evaluate!(y, F.g, x)
    @assert size(y, 1) == size(y, 2)
    @inbounds for i in axes(y, 1)
        y[i, i] += one(eltype(y))
    end
    ldiv!(x, y)
end

function evaluate_diff!(y, F :: GreensFunction, x)
    # (a / x)' = a' / x - a / x^2
    evaluate_diff!(y, F.g, x)
    z = F.g(x)
    axpy!(-inv(x), z, y)
    ldiv(x, y)
end

(F :: GreensFunction)(x) = (F.g(x) + I) / x

# TODO smarter than this
function green_0(oqs :: OQSystem)
    L = Matrix(Commutator(oqs.H, -1))
    Λ, Ψ = eigen(L)
    R = [Ψ[:, i] * Ψ[:, i]' * Λ[i] for i in eachindex(Λ) if !iszero(Λ[i])]
    return GreensFunction(oqs.dim, PoleInterpolation(filter(x -> !iszero(x), Λ), R))
end

function selfconsistency(oqs :: OQSystem, G :: GreensFunction{<: PoleInterpolation}, ω)
    T = promote_type(ComplexF64, eltype(oqs.H))
    ginv = zeros(T, oqs.dim ^ 2, oqs.dim ^ 2)
    tmp2 = zeros(T, oqs.dim ^ 2, oqs.dim ^ 2)
    for b in oqs.baths
        # loop over poles of Green's function
        R = nothing
        K = nothing
        for (p, r) in zip(poles(G.g), residues(G.g))
            if imag(p) > 0.0
                continue
            end
             
            if R isa Nothing
                R = b.R(ω - p)
                K = b.K(ω - p)
            else
                evaluate!(R, b.R, ω - p)
                evaluate!(K, b.K, ω - p)
            end
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(tmp2, left', r, inv(p), zero(T))
                @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                    mul!(ginv, tmp2, right, R[j, k], one(T))
                    mul!(ginv, tmp2, right′, K[j, k], one(T))
                end
            end
        end
        rG = G.g(0.0) + I
        evaluate!(R, b.R, ω)
        evaluate!(K, b.K, ω)
        @inbounds for (j, left) in enumerate(b.cpl_q)
            mul!(tmp2, left', rG)
            @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                mul!(ginv, tmp2, right, R[j, k], one(T))
                mul!(ginv, tmp2, right′, K[j, k], one(T))
            end
        end
    end
    ginv .+= Matrix(Commutator(oqs.H, -1))
    lmul!(-one(T), ginv)
    @inbounds for i in axes(ginv, 1)
        ginv[i, i] += ω
    end
    return inv(ginv)
end

function selfconsistency(oqs :: OQSystem, G :: GreensFunction, ω)
    T = promote_type(ComplexF64, eltype(oqs.H))
    ginv = zeros(T, oqs.dim ^ 2, oqs.dim ^ 2)
    gtmp = Matrix{T}(undef, oqs.dim ^ 2, oqs.dim ^ 2)
    tmp2 = Matrix{T}(undef, oqs.dim ^ 2, oqs.dim ^ 2)
    local rR, rK
    for b in oqs.baths
        # poles of retarded and Keldysh components of
        # bath's Green's function
        @inbounds for (i, p) in enumerate(b.K.poles)
            if i <= length(b.R.poles)
                rR = slice(b.R.residues, i)
            end
            rK = slice(b.K.residues, i)
            evaluate!(gtmp, G, ω - p)
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(tmp2, left', gtmp)
                @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                    if i <= length(b.R.poles)
                        mul!(ginv, tmp2, right, rR[j, k], one(T))
                    end
                    mul!(ginv, tmp2, right′, rK[j, k], one(T))
                end
            end
        end
        # counterterm
        @inbounds for (j, left) in enumerate(b.cpl_q)
            @inbounds for (k, right) in enumerate(b.cpl_c)
                mul!(ginv, left' * right, one(eltype(right)), b.R.cnst[j, k], one(T))
            end
        end
    end
    ginv .+= Matrix(Commutator(oqs.H, -1))
    lmul!(-one(T), ginv)
    @inbounds for i in axes(ginv, 1)
        ginv[i, i] += ω
    end
    return inv(ginv)
end

"""
    steady_state(G :: GreensFunction)

Returns steady state density matrix given converged Green's function.
"""
function steady_state(G :: GreensFunction)
    U, _, _ = svd!(G.g(0.0) + I)
    ret = reshape(U[:, 1], G.dim, G.dim)
    ret ./= tr(ret)
    return ret
end

"""
    simple_iteration(oqs :: OQSystem, G :: GreensFunction, Ω_cutoff :: Real[, which]; kwargs...)

Iterate selfconsistensy equation on Green's function `G` for the system `oqs`.
Returns barycentric or pole interpolation on the segment from `-Ω_cutoff` to `Ω_cutoff`

Keyword arguments:
- `aaa_iter`: number of AAA algorithm iterations, default 20
- `aaa_eps`: accuracy threshold for AAA algorithm, default 1e-9
- `aaa_split`: function which takes iteration index and returns number of splitting points between consecutive support nodes, default `(n -> max(3, 20 - n))`
- `nrm`: norm function, default `LinearAlgebra.norm`
"""
function simple_iteration(
    oqs :: OQSystem,
    G :: GreensFunction,
    Ω_cutoff :: Real;
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)
    F = let oqs = oqs, G = G; ω -> (selfconsistency(oqs, G, ω) * ω - I) end
    F_int = _simple_iteration(oqs, F, Ω_cutoff; aaa_iter, aaa_eps, aaa_split, nrm)

    if G.g isa PoleInterpolation
        return GreensFunction(oqs.dim, PoleInterpolation(F_int))
    end
    return GreensFunction(oqs.dim, F_int)
end