using SelfconsistentBorn
using Test
using LinearAlgebra

function Σ_born(oqs :: OQSystem, ω)
    Λ, Ψ = eigen(oqs.liou)
    iΨ = inv(Ψ)
    return sum(
        sum(
            let R = b.R(ω - λ), K = b.K(ω - λ), l = Ψ[:, i], r = transpose(iΨ[i, :]), cpl_c = b.cpl_c, cpl_q = b.cpl_q
                sum(
                    (left' * l) * (r * (right′ * R[j, k] + right″ * K[j, k]))
                    for (j, left) in enumerate(cpl_q), (k, (right′, right″)) in enumerate(zip(cpl_c, cpl_q))
                ) * (-2π * im)
            end
            for (i, λ) in enumerate(Λ)
        )
        for b in oqs.baths
    ) * im / (2π)
end

@testset "SelfconsistentBorn.jl" begin
    n_fock = 6
    H = Diagonal(0 : n_fock - 1)
    q = SymTridiagonal(zeros(n_fock), [sqrt(i / 2.0) for i in 1 : n_fock - 1])

    bath = BosonicBath([q], ω -> ones(1, 1) * 0.2 * ω / (1.0 + (ω / 20.0) ^ 4), 0.5, 100.0)
    oqs = OQSystem(H, [bath])

    @test norm(Σ_born(oqs, 0.2im) - selfconsistency(oqs, 0.2im)) < 1e-6
    @test norm(Σ_born(oqs, 2.0 + 0.2im) - selfconsistency(oqs, 2.0 + 0.2im)) < 1e-6
    @test norm(Σ_born(oqs, -3.0 + 0.1im) - selfconsistency(oqs, -3.0 + 0.1im)) < 1e-6
end
