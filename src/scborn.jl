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

struct GreensFunction{F <: RationalInterpolation} <: RationalInterpolation
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

(F :: GreensFunction)(x) = (F.g(x) + I) / x

# TODO smarter than this
function green_0(oqs :: OQSystem)
    L = Matrix(Commutator(oqs.H, -1))
    Λ, Ψ = eigen(L)
    return PoleInterpolation(Λ, [Ψ[:, i] * Ψ[:, i]' for i in 1 : size(Ψ, 2)])
end

function green_selfconsistency(oqs :: OQSystem, G, ω)
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
        @inbounds for (k, right) in enumerate(b.cpl_c)
            @inbounds for (j, left) in enumerate(b.cpl_q)
                mul!(ginv, left * right, one(eltype(right)), b.R.cnst[j, k], one(T))
            end
        end
    end
    ginv .+= Matrix(Commutator(oqs.H, -1))
    lmul!(-one(T), ginv)
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
        let oqs = oqs, G = X; ω -> (green_selfconsistency(oqs, G, ω) * ω - I) end
    ) : (
        which == :selfenergy ? nothing : nothing
    )

    # TODO add possibility of a true zero
    ωs = [1e-6, Ω_cutoff]
    fs = reduce(
        append!, F.(ωs);
        init = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, 0)
    )
    ws = [(-1.0 + 0.0im) ^ i for i in 1 : length(ωs)]
    sym_f = conj ∘ oqs.perm
    # symmetric weights and values
    fs′ = similar(fs)
    ws′ = conj(ws)
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
    return GreensFunction(oqs.dim, F_int)
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