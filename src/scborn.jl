struct SCBorn{T, S <: SelfEnergy, B <: BosonicBath}
    H :: T
    Σ :: S
    baths :: Vector{B}
end

green_retarded(s :: SCBorn, ω) = inv(ω * I - s.H - retarded(s.Σ, ω))
green_advanced(s :: SCBorn, ω) = (green_retarded(s, ω'))'
green_advanced(s :: SCBorn, ω :: Real) = (green_retarded(s, ω))'
green_keldysh(s :: SCBorn, ω) = green_retarded(s, ω) * keldysh(s.Σ, ω) * green_advanced(s, ω)
green_keldysh(s :: SCBorn, ω :: Real) = let r = green_retarded(s, ω), k = keldysh(s.Σ, ω)
    r * k * (r')
end

function selfconsistensy(s :: SCBorn{<: Any, SelfEnergy{T, S}, <: Any}, ωs) where {T, S}
    Σ_R = S[]
    Σ_K = S[]
    
    GR = let s = s; ω -> green_retarded(s, ω) end
    #GK = let s = s; ω -> green_keldysh(s, ω) end
    for ω in ωs
        for bath in s.baths
            R(ω′) = let gr = GR(ω′), σk = keldysh(s.Σ, ω′) 
                0.5 * retarded(bath, ω - ω′, gr * σk * (gr')) + keldysh(bath, ω - ω′, gr)
            end
            K(ω′) = let gr = GR(ω′), σk = keldysh(Σ, ω′) 
                0.5 * retarded(bath, ω - ω′, gr) + 
                0.5 * advanced(bath, ω - ω′, gr') +
                keldysh(bath, ω - ω′, gr * σk * (gr'))
            end
            push!(Σ_R, im * sum(
                quadgk(R, a, b)[1] for (a, b) in zip([-Inf; Σ.xs], [Σ.xs; Inf])
            ) / (2π))
            push!(Σ_K, im * sum(
                quadgk(K, a, b)[1] for (a, b) in zip([-Inf; Σ.xs], [Σ.xs; Inf])
            ) / (2π))
        end
    end
    return map(x -> (-0.5im * (x - x')), Σ_R), map(x -> (-0.5im * (x - x')), Σ_K)
end