#= src/snorm.jl

References:

  J. Dixon. Estimating extremal eigenvalues and condition numbers of matrices.
    SIAM J. Numer. Anal. 20 (4): 812-814, 1983.

  J. Kuczynski, H. Wozniakowski. Estimating the largest eigenvalue by the power
    and Lanczos algorithms with a random start. SIAM J. Matrix Anal. Appl. 13
    (4), 1094-1122, 1992.
=#

# spectral norm
function snorm{T}(A::AbstractLinOp{T};
    atol::Real=0, maxiter::Integer=32, rtol::Real=5*eps(real(T)),
    verb::Bool=true)
  atol >= 0 || throw(ArgumentError("atol"))
  rtol >= 0 || throw(ArgumentError("rtol"))
  m, n   = size(A)
  isherm = ishermitian(A)
  xn     = crandn(T, n)
  xm     = Array(T, m)
  xnrm   = vecnorm(xn)
  s      = one(real(T))
  t      = 0
  niter  = 0
  while s > 0 && abs(s - t) > max(atol, t*rtol)
    if niter == maxiter
      verb && warn("iteration limit ($maxiter) reached")
      break
    end
    niter += 1
    scale!(xn, 1/xnrm)
    if isherm
      A_mul_B!(xm, A, xn)
      copy!(xn, xm)
    else
       A_mul_B!(xm, A, xn)
      Ac_mul_B!(xn, A, xm)
    end
    xnrm = vecnorm(xn)
    t = s
    s = isherm ? xnrm : sqrt(xnrm)
  end
  s
end
snorm(A; args...) = snorm(LinOp(A); args...)

# spectral norm difference
snormdiff{T}(A::AbstractLinOp{T}, B::AbstractLinOp{T}; args...) =
  snorm(A - B; args...)
snormdiff(A, B; args...) = snormdiff(LinOp(A), LinOp(B); args...)