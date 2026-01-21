using SelfconsistentBorn
using LinearAlgebra


n_fock = 7
H = Diagonal(0 : n_fock - 1)
q = zeros(n_fock, n_fock)
for i in 1 : n_fock - 1
    q[i, i + 1] = sqrt(i)
end

bath = BosonicBath([(q + q') / sqrt(2.0)], ω -> ones(1, 1) * 0.2 * ω / (1.0 + ω ^ 2 / 400.0), 0.5, 100.0, true; coth = :aaa)
oqs = OQSystem(H, [bath])
Σ = self_energy_0(oqs)
for i in 1 : 5
    global Σ
    Σ = simple_iteration(oqs, Σ, 100.0; aaa_iter = 40)
    display(steady_state(oqs, Σ))
end
