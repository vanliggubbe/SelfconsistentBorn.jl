struct BosonicBath{T, F, G}
    qs :: Vector{T}
    R :: F
    K :: G
end

retarded(b :: BosonicBath, ω) = b.R(ω)
advanced(b :: BosonicBath, ω) = (b.R(ω'))'
advanced(b :: BosonicBath, ω :: Real) = (b.R(ω))'
keldysh(b :: BosonicBath, ω) = b.K(ω)
keldysh(b :: BosonicBath{<: Any, <: Any, <: Real}, ω) = (retarded(b, ω) - advanced(b, ω)) * coth(ω / (2.0 * b.K)) * 0.5
keldysh(b :: BosonicBath{<: Any, <: Any, <: Real}, ω :: Real) = iszero(ω) ? (
    (derivative(b.R, 0.0) |> (x -> (x - x') * 0.5)) * 2.0 * b.K
) : (
    let r = retarded(b, ω), x = coth(ω / (2.0 * b.K))
        (r - r') * x * 0.5
    end
)

couplings(b :: BosonicBath) = b.qs
coupling(b :: BosonicBath, i) = b.qs[i]

for component in (:retarded, :advanced, :keldysh)
    @eval $(component)(b :: BosonicBath, ω, g) = let D = $(component)(b, ω), qs = b.qs
        sum(
            q * g * q′ * D[i, j]
            for (i, q) in enumerate(qs), (j, q′) in enumerate(qs)
        )
    end
end