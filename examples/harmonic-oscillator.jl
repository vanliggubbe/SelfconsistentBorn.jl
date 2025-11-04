using SelfconsistentBorn
using LinearAlgebra
using ProgressBars
using QuadGK

let
    n_fock = 4
    H = Diagonal(0 : n_fock - 1)
    q = zeros(n_fock, n_fock)
    for i in 1 : n_fock - 1
        q[i, i + 1] = sqrt(i)
    end

    DR(ω) = (isfinite(ω) ? (0.2 * ω / (ω + 10.0im)) : (0.2 + 0.0im)) * ones(1, 1)
    bath = BosonicBath([(q + q') / sqrt(2.0)], DR, 0.3)

    G = GreensFunction(
        vcat(
            LinRange(-10, -1, 4),
            LinRange(-1, 4, 18)[2 : end - 1],
            LinRange(4, 10, 4)
        ),
        n_fock
    )
    oqs = SCBorn(H, G, [bath])

    for i_point in 1 : 10
        for i_iter in ProgressBar(1 : 20)
            simple_iteration!(oqs, 0.3)
        end
        println(update_nodes!(oqs, 4))
    end
    ρ = density_matrix(oqs)
    println("Density matrix: ")
    display(ρ)
    println()
    println("Calculated ⟨x̂²⟩ = ", tr(ρ * (q + q') ^ 2 / 2.0))
    GK(ω) = inv(abs2(ω ^ 2 - 1.0 - first(DR(ω)))) * first(keldysh(bath, ω))
    println("     Exact ⟨x̂²⟩ = ", quadgk(GK, -Inf, Inf)[1] * im / (2π))
end