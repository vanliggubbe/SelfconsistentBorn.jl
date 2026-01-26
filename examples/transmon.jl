using SelfconsistentBorn
using LinearAlgebra
using BenchmarkTools

# parameters of transmon
max_q = 6
E_j = 20
C = 0.6
C_c = 0.01
Z = 0.1

H = Diagonal(-max_q : max_q) .^ 2 / (2 * C) - SymTridiagonal(zeros(2 * max_q + 1), ones(2 * max_q)) * E_j / 2
q = Diagonal(-max_q : max_q)

bath = BosonicBath([q], ω -> ones(1, 1) * (ω * C_c ^ 2 * Z) / ((ω * C_c * C * Z) ^ 2 + (C + C_c) ^ 2), 0.5, 400.0, false; coth = :aaa)
oqs = OQSystem(H, (bath,))
Σ0 = self_energy_0(oqs)

Σ = simple_iteration(oqs, Σ0, 400.0, 1.0; aaa_iter = 40)