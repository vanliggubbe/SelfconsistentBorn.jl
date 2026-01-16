module SelfconsistentBorn

import LinearAlgebra: I, norm, tr, axpy!, axpby!, adjoint!, svd!, Diagonal, I, eigvals, lmul!
import QuadGK: quadgk
import ForwardDiff: derivative
import ElasticArrays: ElasticArray, ElasticMatrix, ElasticVector

export GreensFunction, BosonicBath, SCBorn
export retarded, advanced, keldysh, green_retarded, green_advanced, green_keldysh
export density_matrix
export simple_iteration!, update_nodes!

abstract type RationalInterpolation <: Function end
include("utils.jl")

include("aaa.jl")
include("barycentric.jl")
include("pole.jl")

include("bath.jl")

end