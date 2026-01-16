# TODO rewrite better
function aaa_mat_odd(
    a :: Real, b :: Real, f :: Function;
    ε :: Real = 1e-9, n_iter :: Int = 40, split :: Int = 10, λ :: Real = 0.1, Λ :: Real = 1.0
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
        _, j = findmax(i -> norm(fs[i] - gs[i]), js)
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
