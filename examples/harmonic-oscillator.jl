using SelfconsistentBorn
using LinearAlgebra
using BenchmarkTools

n_fock = 8
H = Diagonal(0 : n_fock - 1)
q = SymTridiagonal(zeros(n_fock), sqrt.((1 : n_fock - 1) * 0.5))

bath = BosonicBath([q], ω -> ones(1, 1) * 0.5 * ω / (1.0 + ω ^ 2 / 400.0), 0.2, 1600.0, true; coth = :aaa)

# println("Benchmarking (anti)commutators")
# display(@benchmark mul!($(zeros(ComplexF64, 144, 144)), $(bath.cpl_q[1]'), $(randn(144, 144))))
# display(@benchmark mul!($(zeros(ComplexF64, 144, 144)), $(randn(144, 144)), $(bath.cpl_c[1])))
# display(@benchmark mul!($(zeros(ComplexF64, 144, 144)), $(randn(144, 144)), $(bath.cpl_q[1])))

# creating bath
oqs = OQSystem(H, (bath,))
Σ0 = self_energy_0(oqs)
let
    σ = simple_iteration_domain(oqs)
    Σ = simple_iteration(oqs, Σ0, 1600.0, σ; aaa_iter = 40)
    for i in 0 : 10
        println("=== Line Im ω = $(σ) ===")
        for j in 1 : 10
            Σ = simple_iteration(oqs, Σ, 1600.0, σ; aaa_iter = 40)
            tmp = norm(selfconsistency(oqs, Σ, im * σ) - Σ(im * σ))
            println("!!! Simple iteration error = $(tmp)")
            if tmp < 5e-8
                break
            end
        end
        tmp = norm(selfconsistency(oqs, Σ, 0.0) - Σ(0.0))
        println("### Error = $(tmp)")
        display(eigvals(steady_state(oqs, Σ)))
        poles = SelfconsistentBorn.PoleInterpolation(Σ.Σ).poles .+ im * Σ.shift
        σ = max(0.0, (σ + maximum(imag(poles[imag(poles) .< σ]))) / 2.0)
    end
end
# println("Self-energy evaluation")
# display(@benchmark SelfconsistentBorn.evaluate!($(zeros(ComplexF64, 144, 144)), $Σ0, 1.0))
# println("Selfconsistency evaluation")
# display(@benchmark SelfconsistentBorn.selfconsistency!($zeros(ComplexF64, 144, 144), $oqs, $Σ0, 1.23))

#Σ = simple_iteration(oqs, Σ0, 400.0, 1.0; aaa_iter = 40)
#Σ = simple_iteration(oqs, Σ, 400.0, 0.1; aaa_iter = 40)
# println("Self-energy evaluation")
# display(@benchmark SelfconsistentBorn.evaluate!($(zeros(ComplexF64, 144, 144)), $Σ, 1.0))
# println("Selfconsistency evaluation")
# display(@benchmark SelfconsistentBorn.selfconsistency!($zeros(ComplexF64, 144, 144), $oqs, $Σ, 1.23))

#f = h5open("test.h5", "w")
#f["self-energy"] = Σ
#close(f)
