# TODO rewrite better
function aaa_mat_odd(
    a :: Real, b :: Real, f :: Function;
    weight :: Function = (_ -> 1.0),
    ε :: Real = 1e-9, n_iter :: Int = 60, split :: Int = 10, λ :: Real = 0.1, Λ :: Real = 1.1
)
    xs = collect(LinRange(a, b, split + 2))
    fs = flatten.(f.(xs))
    gs = deepcopy(fs)
    gs .*= 0.0

    XS = eltype(xs)[]
    FS = eltype(fs)[]
    js = collect(1 : length(xs))
    local WS
    for i in 1 : n_iter
        _, j = findmax(i -> norm(fs[i] - gs[i]) * weight(xs[i]), js)
        jj = js[j]
        deleteat!(js, j)
        push!(XS, xs[jj])
        push!(FS, fs[jj])
        left = xs .< XS[end]

        if sum(left) == 0
            if 0.0 < λ < 1.0
                push!(xs, XS[end] * λ)
                push!(fs, flatten(f(xs[end])))
                push!(js, length(xs))
            end
        else
            y = maximum(xs[left])
            append!(xs, LinRange(y, XS[end], split + 2)[2 : end - 1])
            append!(fs, flatten.(f.(xs[end - split + 1 : end])))
            append!(js, length(xs) - split + 1 : length(xs))
        end

        right = xs .> XS[end]
        if sum(right) == 0
            if 1.0 < Λ
                push!(xs, XS[end] * Λ)
                push!(fs, flatten(f(xs[end])))
                push!(js, length(xs))
            end
        else
            y = minimum(xs[right])
            append!(xs, LinRange(XS[end], y, split + 2)[2 : end - 1])
            append!(fs, flatten.(f.(xs[end - split + 1 : end])))
            append!(js, length(xs) - split + 1 : length(xs))
        end

        C = 1 ./ ((xs[js] .^ 2) * ones(i)' - ones(length(js)) * (XS .^ 2)')
        A = ((fs[js] .* xs[js]) * ones(i)' - ones(length(js)) * permutedims(XS .* FS)) .* C
        B = hcat(reduce(vcat, A; dims = 1, init = eltype(first(fs))[])...)
        _, __, V = svd!(B)
        WS = V[:, end]

        gs = deepcopy(fs)
        gs[js] .= (C * (XS .* WS .* FS)) ./ (xs[js] .* (C * WS))
        err = maximum(norm.(fs - gs))
        if err < ε / 2
            break
        end
    end
    shape = size(f(a))
    if !isempty(shape)
        FS = map(x -> reshape(x, shape), FS)
    end
    return OddBarycentricInterpolation(XS[abs.(WS) .> 0.0], WS[abs.(WS) .> 0.0], FS[abs.(WS) .> 0.0])
end

function aaa_real_axis(
    Ω :: Real,
    f :: Function,
    sym :: Function;
    fun_mut :: Bool = false,
    fun_shape :: Tuple{Vararg{Int}} = (),
    fun_type :: Type{<: Number} = Float64,
    sym_mut :: Bool = true,
    aaa_iter :: Integer = 100,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> 3)
)
    @argcheck Ω > 0.0
    @argcheck aaa_iter > 0
    @argcheck aaa_eps > 0.0

    # initialize support nodes
    zs = [Float64(π) / 2.0]
    fs = fun_mut ? let q = ElasticArray{fun_type}(undef, fun_shape..., 1)
        f(slice(q, 1), Ω)
        q
    end : let q = f(Ω);
        ElasticArray(reshape(q, size(q)..., 1))
    end

    # symmetric point
    fs′ = similar(fs)
    if sym_mut
        sym(slice(fs′, 1), slice(fs, 1))
    else
        slice(fs′, 1) .= slice(fs, 1)
    end

    # weight
    ws = [1.0 + 0.0im]

    # initialize intermediate points
    νs = reduce(
        vcat,
        collect(LinRange(a, b, aaa_split(0) + 2)[2 : end - 1])
        #intermediate(a, b, aaa_split(0))
        for (a, b) in zip([0.0;  zs], [zs; Float64(π)])
    )
    gs = ElasticArray{ComplexF64}(undef, firstdims(fs)..., length(νs))
    @inbounds for i in eachindex(νs)
        if fun_mut
            f(slice(gs, i), Ω * cot(νs[i] / 2.0))
        else
            slice(gs, i) .= f(Ω * cot(νs[i] / 2.0))
        end
    end
    tmp2 = zeros(ComplexF64, firstdims(fs)...)
    
    # preparing auxiliary arrays for tall and skinny QR
    nmax = aaa_iter + length(zs)
    R = zeros(Float64, 4 * nmax, 2 * nmax)                                                                  # right triangular matrix from QR
    tmp = Matrix{Float64}(undef, 2 * prod(firstdims(fs)), 2 * nmax)                                         # each layer
    wspace = QRWYWs(length(R) > length(tmp) ? R : tmp; blocksize = min(2 * nmax, 2 * prod(firstdims(fs))))  # LAPACK workspace for QR factorization

    for it in 1 : aaa_iter
        j = 0
        er = -Inf
        @inbounds for i in eachindex(νs)
            fill!(tmp2, zero(eltype(tmp2)))
            den = 0.0im
            for j in eachindex(zs)
                axpy!(ws[j] / (cis(-νs[i]) - cis(-zs[j])), slice(fs, j), tmp2)
                axpy!(conj(ws[j]) / (cis(-νs[i]) - cis(zs[j])), slice(fs′, j), tmp2)
                den += ws[j] / (cis(-νs[i]) - cis(-zs[j])) + conj(ws[j]) / (cis(-νs[i]) - cis(zs[j]))
            end
            ldiv!(den, tmp2)
            axpy!(-1.0, slice(gs, i), tmp2)
            if norm(tmp2) > er
                er = norm(tmp2)
                j = i
            end
        end
        if er < aaa_eps
            break
        end
        @info "Iteration $(it): added node point $(Ω * cot(νs[j] / 2.0)), error $(er)"

        # add a new support point
        new_z = νs[j]
        push!(zs, new_z)
        append!(fs, slice(gs, j))

        # calculate value of symmetric point
        append!(fs′, similar(slice(fs′, 1)))
        if sym_mut
            sym(slice(fs′, lastdim(fs′)), slice(fs, lastdim(fs)))
        else
            slice(fs′, lastdim(fs′)) .= sym(slice(fs, lastdim(fs)))
        end

        # delete points
        left = maximum(zs[zs .< new_z]; init = 0.0)
        right = minimum(zs[zs .> new_z]; init = Float64(π))
        n_left = sum(left .< νs .< new_z)
        n_right = sum(new_z .< νs .< right)
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
            collect(LinRange(left, new_z, max(n_left, aaa_split(it)) + 2)[2 : end - 1]);
            collect(LinRange(new_z, right, max(n_right, aaa_split(it)) + 2)[2 : end - 1]);
            #intermediate(left, new_z, max(n_left, aaa_split(it)));
            #intermediate(new_z, right, max(n_right, aaa_split(it)));
        ]
        resize!(gs, (firstdims(gs)..., cur + length(new_νs)))
        idx = (cur + 1) : (cur + length(new_νs))
        @inbounds for (i, ν) in zip(idx, new_νs)
            if fun_mut
                f(slice(gs, i), Ω * cot(ν / 2.0))
            else
                slice(gs, i) .= f(Ω * cot(ν / 2.0))
            end
        end
        resize!(νs, cur + length(new_νs))
        νs[idx] .= new_νs

        # do tall and skinny SVD
        n = length(zs)                          # number of points, half number of columns
        height_tmp = prod(firstdims(fs)) * 2
        height = min(n * 2, height_tmp)         # height of a single layer

        # construct matrix
        y_ij :: ComplexF64 = 0.0
        z_ij :: ComplexF64 = 0.0
        local f_i, f_j, f̄_j, dst
        height_R :: Int = 0                 # height of current partial QR decomposition
        cols = 1 : (2 * n)
        @inbounds for i in eachindex(νs)
            f_i = vec(slice(gs, i))
            @inbounds for j in eachindex(zs)
                f_j = vec(slice(fs, j))
                f̄_j = vec(slice(fs′, j))
                
                y_ij = inv(cis(-νs[i]) - cis(-zs[j]))
                z_ij = inv(cis(-νs[i]) - cis(zs[j]))

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
        #push!(F_int.perm, length(ws))
        #sortperm!(F_int.perm, abs.(ws))
    end
    qs = expm1.(-im * zs)
    return qs, ws, fs, fs′
    
        #return ω_pol, ω_res
    #ωs = Ω * cot.(zs / 2.0)
    #cnst = -2.0 * real(sum(ws ./ expm1.(-im * zs)))
    #ws′ = 2.0im * Ω * ws ./ expm1.(-im * zs) .^ 2
    #return cnst, ωs, ws′, fs, filter(isfinite, eigvals(A, B))
    #return fs, zs, ws, filter(isfinite, eigvals(A, B))
end

function aaa_real_axis(
    :: Type{<: PoleInterpolation},
    Ω :: Real,
    f :: Function,
    sym :: Function = conj;
    fun_mut :: Bool = false,
    fun_shape :: Tuple{Vararg{Int}} = (),
    fun_type :: Type{<: Number} = Float64,
    sym_mut :: Bool = false,
    aaa_iter :: Integer = 100,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> 3),
    res_eps :: Real = 1e-9
)
    qs, ws, fs, fs′ = aaa_real_axis(Ω, f, sym; aaa_iter, aaa_eps, aaa_split, fun_mut, fun_shape, fun_type, sym_mut)

    A = zeros(2 * length(qs) + 1, 2 * length(qs) + 1)
    B = zeros(2 * length(qs) + 1, 2 * length(qs) + 1)
    for i in eachindex(qs, ws)
        A[2 * i - 1, 2 * i - 1] = real(qs[i])
        A[2 * i - 1, 2 * i - 0] = imag(qs[i])
        A[2 * i - 0, 2 * i - 1] = -imag(qs[i])
        A[2 * i - 0, 2 * i - 0] = real(qs[i])
        A[2 * i - 1, end] = 2 * real(ws[i])
        A[2 * i - 0, end] = -2 * imag(ws[i])
        A[end, 2 * i - 1] = 1.0
        B[2 * i - 1, 2 * i - 1] = 1.0
        B[2 * i - 0, 2 * i - 0] = 1.0
    end
    z_pol = filter(isfinite, eigvals(A, B))
    z_res = [
        sum(
            slice(fs, j) * ws[j] / (z - qs[j]) + slice(fs′, j) * conj(ws[j]) / (z - conj(qs[j]))
            for j in eachindex(ws)
        ) / sum(
            -ws[j] / (z - qs[j]) ^ 2 - conj(ws[j]) / (z - conj(qs[j])) ^ 2
            for j in eachindex(ws)
        )
        for z in z_pol
    ]
    ω_pol = -2im * Ω ./ z_pol .- im * Ω
    ω_res = z_res ./ (1 ./ (ω_pol .+ im * Ω) - (ω_pol .- im * Ω) ./ (ω_pol .+ im * Ω) .^ 2)
    idxs = [norm(slice(ω_res, i)) > res_eps * sqrt(prod(firstdims(ω_res))) for i in axes(ω_res, ndims(ω_res))]
    return PoleInterpolation(ω_pol[idxs], slice(ω_res, idxs), zero(slice(ω_res, 1)))
end

#=
function aaa_real_axis(
    :: Type{<: SymmetricBarycentricInterpolation},
    Ω :: Real,
    f :: Function,
    sym :: Function = conj;
    fun_mut :: Bool = false,
    fun_shape :: Tuple{Vararg{Int}} = (),
    fun_type :: Type{<: Number} = Float64,
    sym_mut :: Bool = false,
    aaa_iter :: Integer = 100,
    aaa_eps :: Real = 1e-9,
    aaa_split :: Function = (n -> 3)
)
    qs, ws, fs, fs′ = aaa_real_axis(Ω, f, sym; aaa_iter, aaa_eps, aaa_split, fun_mut, fun_shape, fun_type, sym_mut)
    ωs = real(-2im * Ω ./ qs .- im * Ω)
    ws′ = 
end
=#
