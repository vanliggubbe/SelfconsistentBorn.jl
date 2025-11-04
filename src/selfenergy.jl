export SelfEnergy

mutable struct SelfEnergy{T <: Real, S}
    ωs :: Vector{T}
    R :: Vector{S}
    K :: Vector{S}
    V :: S

    function SelfEnergy(ωs, R, K, V)
        ωs = collect(ωs)
        if length(ωs) < 2
            throw(ArgumentError("Self-energy approximation must have at least two frequency nodes"))
        end
        if !allunique(ωs)
            throw(ArgumentError("Frequency nodes must be unique"))
        end
        p = sortperm(ωs)
        new{eltype(ωs), promote_type(typeof(V), eltype(R), eltype(K))}(ωs[p], collect(R)[p], collect(K)[p], V)
    end
end

"""
    SelfEnergy(ωs, n :: Int, ν :: Real)

Initialize a frequency independent self-energy for `n`-dimensional system. Retarded component 𝚺ᴿ is given by `-im * ν * I` where I is `n`×`n` identity matrix.
"""
SelfEnergy(ωs, n :: Int, ν :: Real) = SelfEnergy(
    ωs,
    [Matrix(-ComplexF64(ν) * I, n, n) for _ in ωs],
    [Matrix(ComplexF64(2 * ν * (2 - n) / n) * I, n, n) for _ in ωs],
    zeros(ComplexF64, n, n)
)

function add_point!(Σ :: SelfEnergy, ω :: Real)
    push!(Σ.ωs, ω)
    j = length(Σ.ωs)
    while j > 1
        if Σ.ωs[j] < Σ.ωs[j - 1]
            Σ.ωs[j], Σ.ωs[j - 1] = Σ.ωs[j - 1], Σ.ωs[j]
        else
            break
        end
        j -= 1
    end
    if j == 1
        pushfirst!(Σ.R, Σ.R[begin])
        pushfirst!(Σ.K, Σ.K[begin])
    elseif j == length(Σ.ωs)
        push!(Σ.R, Σ.R[end])
        push!(Σ.K, Σ.K[end])
    else
        a, y, b = Σ.ωs[(j - 1) : (j + 1)]
        insert!(Σ.R, j, Σ.R[j - 1] * (b - y) / (b - a) + Σ.R[j] * (y - a) / (b - a))
        insert!(Σ.K, j, Σ.K[j - 1] * (b - y) / (b - a) + Σ.K[j] * (y - a) / (b - a))
    end
end

xlogx(x :: Real) = iszero(x) ? x : (x * log(abs(x)))
xlogx(x :: Integer) = iszero(x) ? 0.0 : (x * log(abs(x)))
xlogx(x :: BigInt) = iszero(x) ? BigFloat(0.0) : (x * log(abs(x)))

function retarded(Σ :: SelfEnergy, x :: Real)
    Re = -(xlogx(x - Σ.ωs[2]) - xlogx(x - Σ.ωs[1])) / (Σ.ωs[2] - Σ.ωs[1]) * Σ.R[1]
    for i in 2 : length(Σ.ωs) - 1
        a, b, c = Σ.ωs[i - 1], Σ.ωs[i], Σ.ωs[i + 1]
        Re .-= (
            xlogx(x - c) / (c - b) + xlogx(x - a) / (b - a) -
            xlogx(x - b) * (c - a) / (c - b) / (b - a)
        ) * Σ.R[i]
    end
    Re .-= (xlogx(x - Σ.ωs[end - 1]) - xlogx(x - Σ.ωs[end])) / (Σ.ωs[end] - Σ.ωs[end - 1]) * Σ.R[end]
    Re ./= π
    Re .+= Σ.V
    if x < Σ.ωs[1]
        return Re + im * Σ.R[1]
    end
    for i in 1 : length(Σ.ωs) - 1
        if Σ.ωs[i] <= x < Σ.ωs[i + 1]
            return Re + im * (
                Σ.R[i] * (Σ.ωs[i + 1] - x) +
                Σ.R[i + 1] * (x - Σ.ωs[i])
            ) / (Σ.ωs[i + 1] - Σ.ωs[i])
        end
    end
    return Re + im * Σ.R[end]
end

advanced(Σ :: SelfEnergy, x :: Real) = retarded(Σ, x)'

function keldysh(Σ :: SelfEnergy, ω :: Real)
    j = searchsortedfirst(Σ.ωs, ω)
    if j == 1
        return im * Σ.K[begin]
    elseif j > length(Σ.ωs)
        return im * Σ.K[end]
    else
        return im * (Σ.K[j - 1] * (Σ.ωs[j] - ω) + Σ.K[j] * (ω - Σ.ωs[j - 1])) / (Σ.ωs[j] - Σ.ωs[j - 1])
    end
end