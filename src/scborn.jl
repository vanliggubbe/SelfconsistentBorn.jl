struct OQSystem{T <: AbstractMatrix, B}
    dim :: Int
    H :: T
    baths :: B

    function OQSystem(H, baths)
        if size(H, 1) != size(H, 2)
            throw(DimensionMismatch("Hamiltonian must be square matrix"))
        end
        dim = size(H, 1)
        for b in baths
            if !(b isa BosonicBath)
                throw(ArgumentError("Each element of `bath` must have a type of bosonic bath"))
            end
        end
        new{typeof(H), typeof(baths)}(dim, H, baths)
    end
end

# TODO smarter than this
function green_0(oqs :: OQSystem)
    L = Matrix(Commutator(oqs.H, -1))
    Λ, Ψ = eigen(L)
    return PoleInterpolation(Λ, [Ψ[:, i] * Ψ[:, i]' for i in 1 : size(Ψ, 2)])
end

function green_selfconsistency(oqs :: OQSystem, G, ω)
    T = promote_type(Complex, eltype(oqs.H))
    ginv = zeros(T, oqs.dim ^ 2, oqs.dim ^ 2)
    gtmp = Matrix{T}(undef, oqs.dim ^ 2, oqs.dim ^ 2)
    tmp2 = Matrix{T}(undef, oqs.dim ^ 2, oqs.dim ^ 2)
    local rR, rK
    for b in oqs.baths
        @inbounds for (i, p) in enumerate(b.K.poles)
            if i <= length(b.R.poles)
                rR = slice(b.R.residues, i)
            end
            rK = slice(b.K.residues, i)
            evaluate!(gtmp, G, ω - p)
            @inbounds for (j, left) in enumerate(b.cpl_q)
                @time mul!(tmp2, left', gtmp)
                @time @inbounds for (k, (right, right′)) in enumerate(zip(b.cpl_c, b.cpl_q))
                    if i <= length(b.R.poles)
                        mul!(ginv, tmp2, right, rR[j, k], one(eltype(ginv)))
                    end
                    mul!(ginv, tmp2, right′, rK[j, k], one(eltype(ginv)))
                end
            end
        end
    end
    ginv .+= Matrix(Commutator(oqs.H, -1))
    lmul!(-one(eltype(ginv)), ginv)
    axpy!(ω, I(oqs.dim ^ 2), ginv)
    return inv(ginv)
end

function simple_iteration(
    oqs :: OQSystem,
    X,
    Ω_cutoff :: Real,
    which = :green;
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)
    F = (which == :green) ? (
        let oqs = oqs, G = X; ω -> green_selfconsistency(oqs, G, ω) end
    ) : (
        which == :selfenergy ? nothing : nothing
    )

    ωs = [-Ω_cutoff, Ω_cutoff]
    fs = reduce(
        append!, F.(ωs);
        init = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 0)
    )
    ws = [(-1.0 + 0.0im) ^ i for i in 1 : length(ωs)]
    F_int = BarycentricInterpolation(ωs, ws, fs)
    @assert F_int.nodes === ωs
    @assert F_int.weights === ws
    @assert F_int.values === fs

    νs = reduce(
        vcat,
        collect(LinRange(a, b, aaa_split(0) + 2))[2 : end - 1]
        for (a, b) in zip(@view(ωs[1 : end - 1]), @view(ωs[2 : end]))
    )
    gs = reduce(
        append!, F.(νs);
        init = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 0)
    )
    for it in 1 : aaa_iter
        er, j = findmax(i -> nrm(F_int(νs[i]) - slice(gs, i)), eachindex(νs))
        if er < aaa_eps
            break
        end
        println(νs[j], " ", er)

        # add a new support point
        new_ω = νs[j]
        push!(ωs, new_ω)
        append!(fs, slice(gs, j))

        # delete points
        left = maximum(ωs[ωs .< new_ω])
        right = minimum(ωs[ωs .> new_ω])
        n_left = sum(left .< νs .< new_ω)
        n_right = sum(new_ω .< νs .< right)
        to_delete = reverse(collect(1 : length(νs))[left .< νs .< right])
        cur = length(νs)
        for x in to_delete
            if x != cur
                slice(gs, x) .= slice(gs, cur)
                νs[x] = νs[cur]
            end
            cur -= 1
        end
        resize!(νs, cur)
        resize!(gs, firstdims(gs)..., cur)

        # add new points
        new_νs = [
            collect(LinRange(left, new_ω, max(n_left, aaa_split(it)) + 2))[2 : end - 1];
            collect(LinRange(new_ω, right, max(n_right, aaa_split(it)) + 2))[2 : end - 1]
        ]
        append!(νs, new_νs)
        for ν in new_νs
            append!(gs, F(ν))
        end

        # do tall and skinny SVD
        n = length(ωs)
        R = zeros(ComplexF64, 2 * n, n)
        tmp = Matrix{ComplexF64}(undef, oqs.dim ^ 4, n)
        wspace = QRWYWs(tmp) # LAPACK workspace for multiple QRs
        @inbounds for i in eachindex(νs)
            tmp2 = vec(slice(gs, i))
            # make a block out of 
            @inbounds for j in eachindex(ωs)
                tmp[:, j] .= tmp2
                axpy!(-one(ComplexF64), vec(slice(fs, j)), view(tmp, :, j))
                ldiv!(νs[i] - ωs[j], view(tmp, :, j))
            end
            q = QRCompactWY(LAPACK.geqrt!(wspace, tmp)...)
            if i == 1
                R[1 : n, :] .= q.R
            else
                R[1 + n : end, :] .= q.R
                LAPACK.geqrt!(wspace, R)
                for j in 1 : n - 1
                    fill!(@view(R[j + 1 : n, j]), zero(eltype(R)))
                end
            end
        end
        _, __, V = svd!(@view R[1 : n, :])
        ws .= V[1 : end - 1, end]
        push!(ws, V[end])
        push!(F_int.perm, length(ws))
        sortperm!(F_int.perm, abs.(ws))
    end
    return F_int
end