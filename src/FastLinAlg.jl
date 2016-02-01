#= src/FastLinAlg.jl
=#

module FastLinAlg

importall Base

export

  # linop.jl
  AbstractLinearOperator,
  LinearOperator,
  HermitianLinearOperator,
  AbstractMatOrLinOp,

  # matrixlib.jl
  matrixlib,

  # snorm.jl
  snorm,
  snormdiff

# source files

include("linop.jl")
include("matrixlib.jl")
include("snorm.jl")
include("util.jl")

end  # module