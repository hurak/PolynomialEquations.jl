module PolynomialEquations

using Polynomials
using LinearAlgebra
using JuMP
using Clp
using Printf


include("transformations.jl")
include("matrices.jl")
include("linear.jl")

export scale
export conjreciprocal
export cconj
export dconj
export ltbtmatrix
export sylvestermatrix
export axby0
export axbyc
export axxabb
export axbycd
export axycminl1

end
