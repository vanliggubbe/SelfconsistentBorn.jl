module SelfconsistentBorn

import LinearAlgebra: I
import QuadGK: quadgk

include("selfenergy.jl")
include("bath.jl")
include("scborn.jl")

end
#=
function update!(bath, Σ, H, q, η)
    V, R, K = selfconsistensy(bath, Σ, Σ.xs, H, q)
    res = max(
        norm(Σ.V - V),
        maximum(norm.(Σ.R - R)),
        maximum(norm.(Σ.K - K))
    )
    Σ.V = (1.0 - η) * Σ.V + V * η
    Σ.R .= (1.0 - η) * Σ.R + R * η
    Σ.K .= (1.0 - η) * Σ.K + K * η
    return res
end

error()
=#