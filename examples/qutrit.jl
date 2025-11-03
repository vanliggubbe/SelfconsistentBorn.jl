using SelfconsistentBorn
using LinearAlgebra

H = Diagonal([0.0, 1.0, 2.1])
q = [0.0 1.0 0.0; 1.0 0.0 1.5; 0.0 1.5 0.0]

DR(ω) = (isfinite(ω) ? (ω / (ω + 10.0im)) : (1.0 + 0.0im)) * ones(1, 1)
bath = BosonicBath([q], DR, 0.01)
Σ = SelfEnergy(
    vcat(
        LinRange(-40, -2, 4),
        LinRange(-2, 2, 10)[2 : end - 1],
        LinRange(2, 40, 4)
    ),
    [Matrix(-Complex(1.0)I, 3, 3) for _ in 1 : 16],
    [zeros(ComplexF64, 3, 3) for _ in 1 : 16],
    zeros(ComplexF64, 3, 3)
)
oqs = SCBorn(H, Σ, [bath])