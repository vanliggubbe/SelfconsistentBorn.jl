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
green_lesser(s :: SCBorn, ω :: Real) = let r = green_retarded(s, ω), k = keldysh(s.Σ, ω)
    (r * k * (r') + (r' - r)) / 2.0
end
green_keldysh_dμ(s :: SCBorn, ω :: Real) = let r = green_retarded(s, ω), k = keldysh(s.Σ, ω)
    -r * (r * k + k * r') * r'
end

function selfconsistensy(s :: SCBorn{<: Any, SelfEnergy{T, S}, <: Any}, ωs; quadgk_kwargs...) where {T, S}
    Σ_R = S[]
    Σ_K = S[]
    
    GR = let s = s; ω -> green_retarded(s, ω) end
    GK = let s = s; ω -> green_keldysh(s, ω) end
    for ω in ωs
        tmp_R = zero(s.Σ.V)
        tmp_K = zero(s.Σ.V)
        for bath in s.baths
            R(ω′) = let gr = GR(ω′), σk = keldysh(s.Σ, ω′) 
                0.5 * retarded(bath, ω - ω′, gr * σk * (gr')) + keldysh(bath, ω - ω′, gr)
            end
            K(ω′) = let gr = GR(ω′), σk = keldysh(s.Σ, ω′) 
                0.5 * retarded(bath, ω - ω′, gr) + 
                0.5 * advanced(bath, ω - ω′, gr') +
                keldysh(bath, ω - ω′, gr * σk * (gr'))
            end
            tmp_R += im * sum(
                quadgk(R, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.Σ.ωs], [s.Σ.ωs; Inf])
            ) / (2π)
            tmp_K += im * sum(
                quadgk(K, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.Σ.ωs], [s.Σ.ωs; Inf])
            ) / (2π)
        end
        push!(Σ_R, tmp_R)
        push!(Σ_K, tmp_K)
    end
    GK_int = im * sum(
        quadgk(GK, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.Σ.ωs], [s.Σ.ωs; Inf])
    ) / (2π)
    return (
        sum(retarded(bath, Inf, GK_int) * 0.5 for bath in s.baths),
        map(x -> (-0.5im * (x - x')), Σ_R),
        map(x -> (-0.5im * (x - x')), Σ_K)
    )
end

density_matrix(s :: SCBorn; quadgk_kwargs...) = -im * sum(
    quadgk(
        let s = s;
            ω -> green_lesser(s, ω)
        end,
        a, b; quadgk_kwargs...
    )[1] for (a, b) in zip([-Inf; s.Σ.ωs], [s.Σ.ωs; Inf])
) / (2π)

density_matrix_dμ(s :: SCBorn; quadgk_kwargs...) = -im * sum(
    quadgk(
        let s = s;
            ω -> green_keldysh_dμ(s, ω)
        end,
        a, b; quadgk_kwargs...
    )[1] for (a, b) in zip([-Inf; s.Σ.ωs], [s.Σ.ωs; Inf])
) / (4π)


function simple_iteration!(s :: SCBorn, η :: Real = 0.5; quadgk_kwargs...)
    if !(0.0 < η <= 1.0)
        throw(ArgumentError("Coefficient in simple iteration must be in (0, 1]"))
    end
    V, R, K = selfconsistensy(s, s.Σ.ωs; quadgk_kwargs...)
    res = max(
        norm(V - s.Σ.V),
        maximum(norm.(s.Σ.R - R)),
        maximum(norm.(s.Σ.K - K))
    )
    s.Σ.V .= (1.0 - η) * s.Σ.V + η * V
    s.Σ.R .= (1.0 - η) * s.Σ.R + η * R
    s.Σ.K .= (1.0 - η) * s.Σ.K + η * K

    tr_ρ = tr(density_matrix(s))
    s.Σ.K .+= s.Σ.R * 4.0 * (tr_ρ - 1.0) / size(s.Σ.V)[1]
    #s.Σ.V -= I * (1.0 - tr_ρ) / (dtr_ρ_dμ)
    return res
end

function update_nodes!(s :: SCBorn, n_split :: Int; quadgk_kwargs...)
    ωs = reduce(
        vcat,
        collect(LinRange(a, b, n_split + 2)[2 : end - 1])
        for (a, b) in zip(@view(s.Σ.ωs[1 : end - 1]), @view(s.Σ.ωs[2 : end]))
    )
    _, R, K = selfconsistensy(s, ωs; quadgk_kwargs...)
    R′ = map(x -> (-0.5im * (x - x')), retarded(s.Σ, ω) for ω in ωs)
    K′ = [-im * keldysh(s.Σ, ω) for ω in ωs]
    err_R = norm.(R - R′)
    err_K = norm.(K - K′)
    j_R = argmax(err_R)
    j_K = argmax(err_K)
    if err_R[j_R] > err_K[j_K]
        add_point!(s.Σ, ωs[j_R])
        return ωs[j_R], err_R[j_R]
    else
        add_point!(s.Σ, ωs[j_K])
        return ωs[j_K], err_K[j_K]
    end
end