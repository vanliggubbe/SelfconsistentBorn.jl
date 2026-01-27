struct OQSystem{T <: AbstractMatrix, B, P}
    dim :: Int
    H :: T
    baths :: B
    perm :: PermutationMap{P}

    L :: Matrix{ComplexF64}
    # allocations for calculations
    ginv :: Matrix{ComplexF64}
    tmp1 :: Matrix{ComplexF64}
    tmp2 :: Matrix{ComplexF64}

    function OQSystem(H, baths)
        if size(H, 1) != size(H, 2)
            throw(DimensionMismatch("Hamiltonian must be square matrix"))
        end
        dim = size(H, 1)
        perm = collect(vec(transpose(reshape(1 : dim ^ 2, dim, dim))))
        for b in baths
            if !(b isa BosonicBath)
                throw(ArgumentError("Each element of `bath` must have a type of bosonic bath"))
            end
            for (c, q) in zip(b.cpl_c, b.cpl_q)
                if size(c) != (dim ^ 2, dim ^ 2) || size(q) != (dim ^ 2, dim ^ 2)
                    throw(DimensionMismatch("Dimensions of the Hamiltonian and the coupling operators must coincide"))
                end
            end
        end
        # liouvillian
        L = Matrix{ComplexF64}(undef, dim ^ 2, dim ^ 2)
        mul!(L, Commutator(H), 1.0)
        new{typeof(H), typeof(baths), typeof(perm)}(dim, H, baths, PermutationMap(perm), L, similar(L), similar(L), similar(L))
    end
end

function _simple_iteration(
    oqs :: OQSystem,
    F :: Function,
    Ω_cutoff :: Real;
    ε :: Real = 1e-12,
    aaa_iter :: Int = 20,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> max(3, 20 - n)),
    nrm :: Function = norm
)

    # initialize support points
    ωs = [Ω_cutoff]
    fs = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, length(ωs))
    @inbounds for i in eachindex(ωs)
        F(slice(fs, i), ωs[i])
    end
    ws = [(-1.0 + 0.0im) ^ i for i in 1 : length(ωs)]

    # symmetric values
    # copy, because permutation is changed later when blocks are identified
    sym_f = conj ∘ PermutationMap(copy(oqs.perm.p))
    fs′ = similar(fs)
    @inbounds for i in axes(fs′, 3)
        evaluate!(slice(fs′, i), sym_f, slice(fs, i))
    end
    F_int = SymmetricBarycentricInterpolation(ωs, ws, fs, conj, sym_f)

    @assert F_int.nodes === ωs
    @assert F_int.weights === ws
    @assert F_int.values === fs

    # points between support
    νs = reduce(
        vcat,
        collect(LinRange(a, b, aaa_split(0) + 2))[2 : end - 1]
        for (a, b) in zip([0.0;  ωs[1 : end - 1]], ωs)
    )
    gs = ElasticArray{ComplexF64}(undef, oqs.dim ^ 2, oqs.dim ^ 2, length(νs))
    @inbounds for i in eachindex(νs)
        F(slice(gs, i), νs[i])
    end

    # preparing auxiliary arrays for tall and skinny QR
    nmax = aaa_iter + length(ωs)
    R = zeros(Float64, 4 * nmax, 2 * nmax)
    tmp = Matrix{Float64}(undef, 2 * oqs.dim ^ 4, 2 * nmax)
    wspace = QRWYWs(length(R) > length(tmp) ? R : tmp; blocksize = min(2 * nmax, 2 * oqs.dim ^ 4)) # LAPACK workspace for multiple QRs
    tmp2 = zeros(ComplexF64, oqs.dim ^ 2, oqs.dim ^ 2)

    nzmask = BitVector(0 for i in 1 : oqs.dim ^ 4)
    for fx in (fs, fs′, gs)
        @inbounds for i in axes(fx, 3)
            @inbounds for j in eachindex(slice(fx, i))
                nzmask[j] |= (abs(slice(fx, i)[j]) > ε)
            end
        end
    end


    for it in 1 : aaa_iter
        # find the point of the grid
        # which is approximated the worst
        j = 0
        er = -Inf
        @inbounds for i in eachindex(νs)
            evaluate!(tmp2, F_int, νs[i])
            axpy!(-1.0, slice(gs, i), tmp2)
            if nrm(tmp2) > er
                er = nrm(tmp2)
                j = i
            end
        end
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
        @inbounds for j in eachindex(slice(fs′, size(fs′, 3)))
            nzmask[j] |= (abs(slice(fs′, size(fs′, 3))[j]) > ε)
        end

        # delete points
        left = maximum(ωs[ωs .< new_ω]; init = 0.0)
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

        # add new points
        new_νs = [
            collect(LinRange(left, new_ω, max(n_left, aaa_split(it)) + 2))[2 : end - 1];
            collect(LinRange(new_ω, right, max(n_right, aaa_split(it)) + 2))[2 : end - 1]
        ]
        resize!(gs, (firstdims(gs)..., cur + length(new_νs)))
        idx = (cur + 1) : (cur + length(new_νs))
        @inbounds for (i, ν) in zip(idx, new_νs)
            F(slice(gs, i), ν)
            @inbounds for j in eachindex(slice(gs, i))
                nzmask[j] |= (abs(slice(gs, i)[j]) > ε)
            end
        end
        resize!(νs, cur + length(new_νs))
        νs[idx] .= new_νs

        # do tall and skinny SVD
        n = length(ωs)                      # number of points, half number of columns
        idxs = (1 : oqs.dim ^ 4)[nzmask]    # pick only elements with non-zero values
        height_tmp = length(idxs) * 2
        height = min(n, length(idxs)) * 2   # height of a single layer
        
        y_ij :: ComplexF64 = 0.0
        z_ij :: ComplexF64 = 0.0
        local f_i, f_j, f̄_j, dst
        height_R :: Int = 0                 # height of current partial QR decomposition
        cols = 1 : (2 * n)
        @inbounds for i in eachindex(νs)
            f_i = view(vec(slice(gs, i)), idxs)
            @inbounds for j in eachindex(ωs)
                f_j = view(vec(slice(fs, j)), idxs)
                f̄_j = view(vec(slice(fs′, j)), idxs)
                
                y_ij = inv(νs[i] - ωs[j])
                z_ij = inv(νs[i] + ωs[j])

                # real values
                dst = reinterpret(ComplexF64, @view tmp[1 : height_tmp, 2 * j - 1])
                mul!(dst, y_ij, f_j)
                axpy!(z_ij, f̄_j, dst)
                axpy!(-y_ij - z_ij, f_i, dst)
                # imag values
                dst = reinterpret(ComplexF64, @view tmp[1 : height_tmp, 2 * j])
                mul!(dst, im * y_ij, f_j)
                axpy!(-im * z_ij, f̄_j, dst)
                axpy!(-im * y_ij + im * z_ij, f_i, dst)
            end
            # do QR factorization of the slice
            LAPACK.geqrt!(wspace, view(tmp, 1 : height_tmp, cols))
            triu!(@view tmp[1 : height, cols])
            if i == 1
                @views R[1 : height, cols] .= tmp[1 : height, cols]
                height_R = height
            else
                @views R[(height_R + 1) : (height_R + height), cols] .= tmp[1 : height, cols]
                height_R += height
                LAPACK.geqrt!(wspace, @view R[1 : height_R, cols])
                height_R = min(height_R, 2 * n)
                triu!(@view R[1 : height_R, cols])
            end
        end
        _, __, V = svd!(@view R[1 : height_R, cols])
        ws .= reinterpret(ComplexF64, @view V[1 : end - 2, end])
        push!(ws, V[end - 1] + im * V[end])
        push!(F_int.perm, length(ws))
        sortperm!(F_int.perm, abs.(ws))
    end
    return F_int
end

include("greens.jl")
include("selfenergy.jl")