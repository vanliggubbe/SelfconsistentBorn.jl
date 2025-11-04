mutable struct GreensFunction{T <: Real, S}
    ωs :: Vector{T}
    R :: Vector{S}
    K :: Vector{S}

    function GreensFunction(ωs, R, K)
        ωs = sort(collect(ωs))
        if length(ωs) < 3
            throw(ArgumentError("Self-energy approximation must have at least two frequency nodes"))
        end
        if !allunique(ωs)
            throw(ArgumentError("Frequency nodes must be unique"))
        end
        if !allequal([length(ωs), length(R) + 2, length(K) + 2])
            throw(ArgumentError("Length mismatch"))
        end
        new{eltype(ωs), promote_type(eltype(R), eltype(K))}(ωs, R, K)
    end
end

_R(G :: GreensFunction, i) = (i <= 1 || i >= length(G.ωs)) ? zero(G.R[begin]) : G.R[i - 1]
_K(G :: GreensFunction, i) = (i <= 1 || i >= length(G.ωs)) ? zero(G.K[begin]) : G.K[i - 1]

"""
    GreensFunction(ωs)

Make a dummy Green's function with frequency nodes given by `ωs`
"""
GreensFunction(ωs, n :: Int) = normalize!(GreensFunction(
    ωs, [zeros(ComplexF64, n, n) for _ in 1 : length(ωs) - 2],
    [zeros(ComplexF64, n, n) for _ in 1 : length(ωs) - 2]
))

function normalize!(G :: GreensFunction)
    K = tr(int_K(G))
    RA = -2π * im * I - int_RA(G)
    sω = im * (-G.ωs[2] - G.ωs[1] + G.ωs[end] + G.ωs[end - 1])
    D = size(RA)[1]
    for i in 1 : length(G.ωs) - 2
        G.R[i] .+= RA / sω
        G.K[i] .-= (1 .- D / 2.0 + im * K / (4π)) / (im * D * sω / (8π)) * I(D)
    end
    return G
end

function add_point!(G :: GreensFunction, ω :: Real)
    push!(G.ωs, ω)
    j = length(G.ωs)
    while j > 1
        if G.ωs[j] < G.ωs[j - 1]
            G.ωs[j], G.ωs[j - 1] = G.ωs[j - 1], G.ωs[j]
        else
            break
        end
        j -= 1
    end
    if j == 1
        pushfirst!(G.R, zero(G.R[begin]))
        pushfirst!(G.K, zero(G.K[begin]))
    elseif j == length(G.ωs)
        push!(G.R, zero(G.R[end]))
        push!(G.K, zero(G.K[end]))
    else
        a, y, b = G.ωs[(j - 1) : (j + 1)]
        insert!(G.R, j - 1, zero(G.R[begin]))
        insert!(G.K, j - 1, zero(G.R[begin]))
        if j - 2 >= 1
            G.R[j - 1] .+= G.R[j - 2] * ((b - y) / (b - a))
            G.K[j - 1] .+= G.K[j - 2] * ((b - y) / (b - a))
        end
        if j - 1 <= length(G.R)
            G.R[j - 1] .+= G.R[j - 1] * ((y - a) / (b - a))
            G.K[j - 1] .+= G.K[j - 1] * ((y - a) / (b - a))
        end
    end
end

xlogx(x :: Real) = iszero(x) ? x : (x * log(abs(x)))
xlogx(x :: Integer) = iszero(x) ? 0.0 : (x * log(abs(x)))
xlogx(x :: BigInt) = iszero(x) ? BigFloat(0.0) : (x * log(abs(x)))

function retarded(G :: GreensFunction, x :: Real)
    Re = zero(G.R[begin])
    for i in 2 : length(G.ωs) - 1
        a, b, c = G.ωs[i - 1], G.ωs[i], G.ωs[i + 1]
        Re .-= (
            xlogx(x - c) / (c - b) + xlogx(x - a) / (b - a) -
            xlogx(x - b) * (c - a) / (c - b) / (b - a)
        ) * G.R[i - 1]
    end
    Re ./= π
    if x < G.ωs[1]
        return Re
    end
    for i in 1 : length(G.ωs) - 1
        if G.ωs[i] <= x < G.ωs[i + 1]
            return Re + im * (
                _R(G, i) * (G.ωs[i + 1] - x) +
                _R(G, i + 1) * (x - G.ωs[i])
            ) / (G.ωs[i + 1] - G.ωs[i])
        end
    end
    return Re
end

advanced(G :: GreensFunction, x :: Real) = retarded(G, x)'

function keldysh(G :: GreensFunction, ω :: Real)
    j = searchsortedfirst(G.ωs, ω)
    if j == 1
        return im * zero(G.K[begin])
    elseif j > length(G.ωs)
        return im * zero(G.K[end])
    else
        return im * (_K(G, j - 1) * (G.ωs[j] - ω) + _K(G, j) * (ω - G.ωs[j - 1])) / (G.ωs[j] - G.ωs[j - 1])
    end
end

int_K(G :: GreensFunction) = 0.5im * sum(g * (G.ωs[i + 2] - G.ωs[i]) for (i, g) in enumerate(G.K))
int_RA(G :: GreensFunction) = im * sum(g * (G.ωs[i + 2] - G.ωs[i]) for (i, g) in enumerate(G.R))