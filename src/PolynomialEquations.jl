module PolynomialEquations

using Polynomials
using LinearAlgebra

include("transformations.jl")
include("matrices.jl")
include("linear.jl")

export scale
export reciprocal
export paraconj
export ltbtmatrix
export sylvestermatrix
export axbyc
export axby0

end
