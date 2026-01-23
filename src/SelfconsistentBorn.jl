module SelfconsistentBorn

import LinearAlgebra: I, norm, tr, axpy!, axpby!, adjoint!, svd!, Diagonal, eigvals, lmul!, mul!, eigen, ldiv!, QRCompactWY, tr, LU, rdiv!, svd, ishermitian, issymmetric
import ElasticArrays: ElasticArray, ElasticMatrix, ElasticVector
import LinearMaps: LinearMap, _unsafe_mul!, MulStyle, FiveArg, AdjointMap, TransposeMap
import FastLapackInterface: QRWYWs, LAPACK, LUWs
import SpecialFunctions: digamma

export BosonicBath, OQSystem
export green_0, self_energy_0, selfconsistency, simple_iteration, steady_state
export Commutator

abstract type RationalInterpolation <: Function end
include("utils.jl")

include("aaa.jl")
include("barycentric.jl")
include("pole.jl")
include("commutators.jl")
include("permutations.jl")

include("bath.jl")
include("scborn.jl")

end