using SelfconsistentBorn
using LinearAlgebra


n_fock = 4
H = Diagonal(0 : n_fock - 1)
q = zeros(n_fock, n_fock)
for i in 1 : n_fock - 1
    q[i, i + 1] = sqrt(i)
end

bath = BosonicBath([(q + q') / sqrt(2.0)], ω -> ones(1, 1) * 0.1 * ω / (1.0 + ω ^ 2 / 100.0), 0.1, 30.0)
oqs = OQSystem(H, [bath])
g0 = green_0(oqs)
