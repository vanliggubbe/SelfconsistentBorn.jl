using SelfconsistentBorn
using LinearAlgebra
using LinearMaps

n_fock = 12
H = Diagonal(0 : n_fock - 1)
q = SymTridiagonal(zeros(n_fock), sqrt.((1 : n_fock - 1) * 0.5))

bath = BosonicBath([q], ω -> ones(1, 1) * 0.05 * ω / (1.0 + ω ^ 2 / 400.0), 0.8, 400.0, true; coth = :aaa)
oqs = OQSystem(H, [bath])
Σ = self_energy_0(oqs)
for i in 1 : 1
    global Σ
    Σ = simple_iteration(oqs, Σ, 400.0, 6.0 - i; aaa_iter = 40)
    display(steady_state(oqs, Σ))
end
