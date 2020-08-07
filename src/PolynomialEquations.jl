module PolynomialEquations

using Polynomials
using LinearAlgebra
using JuMP
using Clp
using Printf


include("transformations.jl")
include("matrices.jl")
include("linear.jl")
include("analysis.jl")

export scale
export conjreciprocal
export cconj
export dconj
export ltbtmatrix
export sylvesterresultantmatrix
export sylvestermatrix
export iscoprime
export isschurstable
export ishurwitzstable
export axb
export axby0
export axbyc
export axaxbb
export axbycd
export axycminl1

end
