using SelfconsistentBorn
using LinearAlgebra


n_fock = 6
H = Diagonal(0 : n_fock - 1)
q = zeros(n_fock, n_fock)
for i in 1 : n_fock - 1
    q[i, i + 1] = sqrt(i)
end

bath = BosonicBath([(q + q') / sqrt(2.0)], ω -> ones(1, 1) * 0.05 * ω / (1.0 + ω ^ 2 / 400.0), 0.5, 60.0, true)
oqs = OQSystem(H, [bath])
G = green_0(oqs)
for i in 1 : 5
    global G
    G = simple_iteration(oqs, G, 60.0; aaa_iter = 30)
end
