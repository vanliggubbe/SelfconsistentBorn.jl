struct SCBorn{T, S <: GreensFunction, B <: BosonicBath}
    H :: T
    G :: S
    baths :: Vector{B}
end

green_retarded(s :: SCBorn, ω) = retarded(s.G, ω)
green_advanced(s :: SCBorn, ω) = advanced(s.G, ω)
green_keldysh(s :: SCBorn, ω) = keldysh(s.G, ω)
green_lesser(s :: SCBorn, ω :: Real) = let r = green_retarded(s, ω), k = green_keldysh(s, ω)
    (r - r' - k) / 2.0
end

function selfconsistensy(s :: SCBorn{<: Any, GreensFunction{T, S}, <: Any}, ωs; quadgk_kwargs...) where {T, S}
    G_R = S[]
    G_K = S[]
    
    for ω in ωs
        Σ_R = nothing
        Σ_K = nothing
        for bath in s.baths
            R(ω′) = let r = green_retarded(s, ω′), k = green_keldysh(s, ω′) 
                0.5 * retarded(bath, ω - ω′, k) + keldysh(bath, ω - ω′, r)
            end
            K(ω′) = let r = green_retarded(s, ω′), k = green_keldysh(s, ω′) 
                0.5 * spectral_dos(bath, ω - ω′, r - r') +
                keldysh(bath, ω - ω′, k)
            end
            if Σ_R isa Nothing
                Σ_R = im * sum(
                    quadgk(R, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.G.ωs], [s.G.ωs; Inf])
                ) / (2π)
                Σ_K = im * sum(
                    quadgk(K, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.G.ωs], [s.G.ωs; Inf])
                ) / (2π)
            else
                Σ_R .+= im * sum(
                    quadgk(R, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.G.ωs], [s.G.ωs; Inf])
                ) / (2π)
                Σ_K .+= im * sum(
                    quadgk(K, a, b; quadgk_kwargs...)[1] for (a, b) in zip([-Inf; s.G.ωs], [s.G.ωs; Inf])
                ) / (2π)
            end
            
        end
        #display(Σ_R)
        #display(Σ_K)
        push!(G_R, inv(ω * I - s.H - Σ_R))
        push!(G_K, G_R[end] * Σ_K * G_R[end]')
    end
    #display(G_R[1])
    #display(G_K[1])
    return (
        map(x -> (-0.5im * (x - x')), G_R),
        map(x -> (-0.5im * (x - x')), G_K)
    )
end

density_matrix(s :: SCBorn) = im * (int_RA(s.G) - int_K(s.G)) / (4π)

function simple_iteration!(s :: SCBorn, η :: Real = 0.5; quadgk_kwargs...)
    if !(0.0 < η <= 1.0)
        throw(ArgumentError("Coefficient in simple iteration must be in (0, 1]"))
    end
    R, K = selfconsistensy(s, s.G.ωs[2 : end - 1]; quadgk_kwargs...)
    res = max(
        maximum(norm.(s.G.R - R)),
        maximum(norm.(s.G.K - K))
    )
    s.G.R .= (1.0 - η) * s.G.R + η * R
    s.G.K .= (1.0 - η) * s.G.K + η * K

    normalize!(s.G)
    return res
end

function update_nodes!(s :: SCBorn, n_split :: Int; quadgk_kwargs...)
    ωs = reduce(
        vcat,
        collect(LinRange(a, b, n_split + 2)[2 : end - 1])
        for (a, b) in zip(@view(s.G.ωs[1 : end - 1]), @view(s.G.ωs[2 : end]))
    )
    R, K = selfconsistensy(s, ωs; quadgk_kwargs...)
    R′ = map(x -> (-0.5im * (x - x')), retarded(s.G, ω) for ω in ωs)
    K′ = [-im * keldysh(s.G, ω) for ω in ωs]
    err_R = norm.(R - R′)
    err_K = norm.(K - K′)
    j_R = argmax(err_R)
    j_K = argmax(err_K)
    if err_R[j_R] > err_K[j_K]
        add_point!(s.G, ωs[j_R])
        return ωs[j_R], err_R[j_R]
    else
        add_point!(s.G, ωs[j_K])
        return ωs[j_K], err_K[j_K]
    end
end