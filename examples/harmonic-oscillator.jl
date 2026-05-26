using SelfconsistentBorn
using LinearAlgebra
using BenchmarkTools

n_fock = 6
H = Diagonal(0 : n_fock - 1)
q = SymTridiagonal(zeros(n_fock), sqrt.((1 : n_fock - 1) * 0.5))

bath = BosonicBath([q], ω -> ones(1, 1) * 0.01 * ω / (1.0 + ω ^ 2 / 400.0) / (1.0 + ω ^ 2 / 1600), 1.0, 100.0, true)

# creating bath
oqs = OQSystem(H, [bath,])
simple_iteration!(oqs, 100; simple_iter = 10)
σ = simple_iteration_domain(oqs)
for i in 1 : 20
    global σ
    simple_iteration!(oqs, 100, σ; simple_iter = 5)
    println("σ = $(σ), selfconsistency check:", selfconsistency_check(oqs, σ * im), " ", selfconsistency_check(oqs, 0.0))
    display(eigvals(steady_state(oqs)))
    σ *= 0.8
end

