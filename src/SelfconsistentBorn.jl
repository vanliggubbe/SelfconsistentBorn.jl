module SelfconsistentBorn

import LinearAlgebra: I, norm, tr, axpy!, axpby!, adjoint!, svd!, Diagonal, eigvals, lmul!, mul!, eigen, ldiv!, QRCompactWY
import ElasticArrays: ElasticArray, ElasticMatrix, ElasticVector
import LinearMaps: LinearMap, _unsafe_mul!, Adjoint
import FastLapackInterface: QRWYWs, LAPACK

export BosonicBath, OQSystem
export green_0, green_selfconsistency, simple_iteration

abstract type RationalInterpolation <: Function end
include("utils.jl")

include("aaa.jl")
include("barycentric.jl")
include("pole.jl")
include("commutators.jl")

include("bath.jl")
include("scborn.jl")

end