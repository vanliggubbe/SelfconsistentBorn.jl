using SelfconsistentBorn
using Test
using LinearAlgebra
using HDF5
using SparseArrays

@testset "SelfconsistentBorn.jl" begin
    n_fock = 6
    H = Diagonal(0 : n_fock - 1)
    q = SymTridiagonal(zeros(n_fock), [sqrt(i / 2.0) for i in 1 : n_fock - 1])

    bath = BosonicBath([q], ω -> ones(1, 1) * 0.05 * ω / (1.0 + ω ^ 2 / 400.0), 0.5, 100.0, true; coth = :aaa)
    oqs = OQSystem(H, [bath])
    Σ = self_energy_0(oqs)
    Σ = simple_iteration(oqs, Σ, 100.0, 1.0; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.5; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.25; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.125; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.125; aaa_iter = 30)


    for ω in -5 : 0.5 : 5
        @test norm(Σ(ω) - selfconsistency(oqs, Σ, ω)) < 1e-2
    end

    @testset "IO with HDF5" begin
        file_as_bytes = h5open("self_energy", "w"; driver = Drivers.Core(; backing_store = false)) do f
            f["self-energy"] = Σ
            return Vector{UInt8}(f)
        end

        f = h5open(file_as_bytes, "r"; name = "self_energy.h5")
        Σ_new = read(f, "self-energy", SelfconsistentBorn.SelfEnergy)
        close(f)
        for ω in -5 : 0.5 : 5
            @test Σ(ω) ≈ Σ_new(ω)
        end
    end

    max_q = 5
    q1 = sparse(circshift(Matrix(1.0I, 2 * max_q, 2 * max_q), (0, 1)))
    q2 = sparse(circshift(Matrix(1.0I, 2 * max_q, 2 * max_q), (0, -1)))
    q3 = Diagonal(-max_q : max_q - 1)
    bath1 = BosonicBath((q1, q2), ω -> Matrix(1.0I, 2, 2) * 0.05 * ω / (1.0 + ω ^ 2 / 400.0), 0.5, 100.0, false; coth = :aaa)
    bath2 = BosonicBath((q3,), ω -> ones(1, 1) * 0.01 * ω / (1.0 + ω ^ 2 / 100.0), 1.0, 100.0, false; coth = :aaa)

    H = Diagonal((-max_q : max_q - 1) .^ 2) - 2.0 * (q1 * q1 + q2 * q2)
    oqs = OQSystem(H, (bath1, bath2))

    Σ = self_energy_0(oqs)
    Σ = simple_iteration(oqs, Σ, 100.0, 1.0; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.5; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.25; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.125; aaa_iter = 30)
    Σ = simple_iteration(oqs, Σ, 100.0, 0.125; aaa_iter = 30)

    for ω in -5 : 0.5 : 5
        @test norm(Σ(ω) - selfconsistency(oqs, Σ, ω)) < 1e-2
    end

    @testset "IO with HDF5" begin
        file_as_bytes = h5open("self_energy", "w"; driver = Drivers.Core(; backing_store = false)) do f
            f["self-energy"] = Σ
            return Vector{UInt8}(f)
        end

        f = h5open(file_as_bytes, "r"; name = "self_energy.h5")
        Σ_new = read(f, "self-energy", SelfconsistentBorn.SelfEnergy)
        close(f)
        for ω in -5 : 0.5 : 5
            @test Σ(ω) ≈ Σ_new(ω)
        end
    end
end
