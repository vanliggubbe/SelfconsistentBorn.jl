struct OQSystem{T <: AbstractMatrix, B, P}
    dim :: Int
    H :: T
    baths :: B
    perm :: PermutationMap{P}

    function OQSystem(H, baths)
        if size(H, 1) != size(H, 2)
            throw(DimensionMismatch("Hamiltonian must be square matrix"))
        end
        dim = size(H, 1)
        perm = vec(transpose(reshape(1 : dim ^ 2, dim, dim)))
        for b in baths
            if !(b isa BosonicBath)
                throw(ArgumentError("Each element of `bath` must have a type of bosonic bath"))
            end
            if firstdims(b.cpl) != (dim, dim)
                throw(DimensionMismatch("Dimensions of the Hamiltonian and the coupling operators must coincide"))
            end
        end
        new{typeof(H), typeof(baths), typeof(perm)}(dim, H, baths, PermutationMap(perm))
    end
end

function _simple_iteration(
    oqs :: OQSystem,
    F :: Function,
    Ω_cutoff :: Real;
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)
    # TODO add possibility of a true zero
    ωs = [1e-3, Ω_cutoff]
    fs = reduce(
        append!, F.(ωs);
        init = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 0)
    )
    ws = [(-1.0 + 0.0im) ^ i for i in 1 : length(ωs)]
    sym_f = conj ∘ oqs.perm
    # symmetric values
    fs′ = similar(fs)
    @inbounds for i in axes(fs′, 3)
        evaluate!(slice(fs′, i), sym_f, slice(fs, i))
    end
    F_int = SymmetricBarycentricInterpolaton(ωs, ws, fs, conj, sym_f)
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
        @info "Iteration $(it): added node point $(νs[j]), error $(er)"

        # add a new support point
        new_ω = νs[j]
        push!(ωs, new_ω)
        append!(fs, slice(gs, j))
        # calculate value of symmetric point
        append!(fs′, similar(slice(fs′, 1)))
        evaluate!(slice(fs′, size(fs′, 3)), sym_f, slice(fs, size(fs, 3)))

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
        R = zeros(Float64, 4 * n, 2 * n)
        tmp = Matrix{Float64}(undef, 2 * oqs.dim ^ 4, 2 * n)
        height = min(n, oqs.dim ^ 4) * 2
        wspace = QRWYWs(tmp) # LAPACK workspace for multiple QRs

        local y_ij, z_ij, f_i, f_j, f̄_j
        @inbounds for i in eachindex(νs)
            f_i = vec(slice(gs, i))
            # views on real and imaginary parts
            #mul!(tmp, -one(Float64), reinterpret(Float64, reshape(fs, oqs.dim ^ 4, n)))
            @inbounds for j in eachindex(ωs)
                f_j = vec(slice(fs, j))
                f̄_j = vec(slice(fs′, j))
                
                y_ij = inv(νs[i] - ωs[j])
                z_ij = inv(νs[i] + ωs[j])

                # real values
                dst = reinterpret(ComplexF64, @view tmp[:, 2 * j - 1])
                mul!(dst, y_ij, f_j)
                axpy!(z_ij, f̄_j, dst)
                axpy!(-y_ij - z_ij, f_i, dst)
                # imag values
                dst = reinterpret(ComplexF64, @view tmp[:, 2 * j])
                mul!(dst, im * y_ij, f_j)
                axpy!(-im * z_ij, f̄_j, dst)
                axpy!(-im * y_ij + im * z_ij, f_i, dst)
            end
            q = QRCompactWY(LAPACK.geqrt!(wspace, tmp)...)
            if i == 1
                R[1 : height, :] .= q.R
            else
                R[n .+ (1 : height), :] .= q.R
                LAPACK.geqrt!(wspace, R)
                for j in 1 : 2 * n
                    fill!(@view(R[j + 1 : end, j]), zero(eltype(R)))
                end
            end
        end
        _, __, V = svd!(@view R[1 : 2 * n, :])
        ws .= reinterpret(ComplexF64, @view V[1 : end - 2, end])
        push!(ws, V[end - 1] + im * V[end])
        push!(F_int.perm, length(ws))
        sortperm!(F_int.perm, abs.(ws))
    end
    return F_int
end

include("greens.jl")
include("selfenergy.jl")