module SelfconsistentBorn

import LinearAlgebra: I, norm, tr, axpy!, axpby!, adjoint!, svd!, Diagonal, eigvals, lmul!, rmul!, mul!, eigen, ldiv!, QRCompactWY, tr, LU, rdiv!, svd, ishermitian, issymmetric, triu!, svdvals
import ElasticArrays: ElasticArray, ElasticMatrix, ElasticVector
import LinearMaps: LinearMap, _unsafe_mul!, MulStyle, FiveArg, AdjointMap, TransposeMap
import FastLapackInterface: QRWYWs, LAPACK, LUWs
import SpecialFunctions: digamma
import DataStructures: IntDisjointSet, num_groups, find_root!
import HDF5: File, Group, create_group, attributes, read_attribute
import ArgCheck: @argcheck
import SimpleTraits: @traitdef, @traitimpl, @traitfn

export BosonicBath, OQSystem
export selfconsistency, simple_iteration!, steady_state, selfconsistency_check, dim, self_energy
export Commutator
export simple_iteration_domain

include("utils.jl")

abstract type RationalInterpolation{N} <: Function end

include("barycentric.jl")
include("pole.jl")
include("aaa.jl")
include("commutators.jl")
include("permutations.jl")

include("bath.jl")
include("scborn.jl")

end
