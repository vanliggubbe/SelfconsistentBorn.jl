using SelfconsistentBorn
using LinearAlgebra

n_fock = 4
H = Diagonal(0 : n_fock - 1)
q = zeros(n_fock, n_fock)
for i in 1 : n_fock - 1
    q[i, i + 1] = sqrt(i)
end

DR(ω) = (isfinite(ω) ? (0.2 * ω / (ω + 10.0im)) : (0.2 + 0.0im)) * ones(1, 1)
bath = BosonicBath([q + q'], DR, 0.5)

Σ = SelfEnergy(
    vcat(
        LinRange(-40, -2, 4),
        LinRange(-2, 2, 10)[2 : end - 1],
        LinRange(2, 40, 4)
    ),
    n_fock, 0.5
)
oqs = SCBorn(H, Σ, [bath])