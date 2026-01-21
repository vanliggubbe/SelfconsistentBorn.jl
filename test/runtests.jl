using SelfconsistentBorn
using Test
using LinearAlgebra

@testset "SelfconsistentBorn.jl" begin
    n_fock = 3
    H = Diagonal(0 : n_fock - 1)
    q = zeros(n_fock, n_fock)
    for i in 1 : n_fock - 1
        q[i, i + 1] = sqrt(i)
    end

    bath = BosonicBath([(q + q') / sqrt(2.0)], ω -> ones(1, 1) * 0.2 * ω / (1.0 + ω ^ 2 / 400.0), 0.5, 100.0, true; coth = :digamma)
    oqs = OQSystem(H, [bath])
    G0 = green_0(oqs)
    G0_brute(ω) = inv(ω * I - Matrix(SelfconsistentBorn.Commutator(H, -1)))
    @test G0(0.5) ≈ G0_brute(0.5)
    @test G0(-1.3) ≈ G0_brute(-1.3)
end
