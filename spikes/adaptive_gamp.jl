# This implements Adaptive Generalized Approximate Message Parsing (adGAMP)
# for compressed sensing / regularized regression as described in:
#   Kamilov et al, "Approximate Message Passing with Consistent Parameter
#   Estimation and Applications to Sparse Learning", arxiv.org, 2012,
#   http://arxiv.org/pdf/1207.3859v3.pdf
#

# Unpack the adGamp parameters and options.
#
# Output:
#  Tuple of parameters
#
# Input:
#  A - measurement matrix (M*N) where M is the number of cases and N the number of covariates
#  y - noisy measurements (M*1) assuming y = A * x + err
#
# Options:
#  x0 - starting point for search (if used in sequential scenarios etc)
#  maxIterations -
#  maxParamIterations -
function unpack_adGamp_params(A, y;
  adapt = true,
  maxIterations      = 100,
  maxParamIterations = 300,
  x0 = nothing,       # A starting point, will be set to zero if not given

  f::Function = scad, # Defining function used in tresholding, by default SCAD. Should be soft?

  tolerance = 1e-6,   # Convergence level, i.e. terminate algorithm if norm change is less than this between iterations

  )

  M, N = size(A)
  xhat = x0 || 1

end

  a = 1

immutable AdGampOptions
  adapt::Bool
  maxIterations::Int64
  maxParamIterations::Int64
end

function adaptiveGamp(A, y; kws...)
  A, y, denoiseFunc, estimatingFuncs, xhat, vx, theta, options = unpack_adGamp_params(A, y; kws...)

  m, n = size(A)
  trA = A'
  hadamardA = A .* A
  trhadamardA = hadamardA'
  xhatprev = 1e6 * ones(n, 1)
  shat = zeros(m, 1)
  numBelowTolerance = 0

  for i in 1:maxIterations
    vp = hadamardA * vx
    phat = A * xhat .- (vp .* shat)
    vs = 1 ./ (vp + noiseVariance)
    shat = vs .* (y .- phat)

    # Update variables
    vr = 1 ./ (trhadamardA * vs)
    rhat = xhat .+ (vr .* (trA * shat))

    if options.adapt
      theta = update_distributional_params(rhat, vr, theta, estimatingFuncs, options)
    end

    xhat, vx = denoiseFunc(rhat, vr, theta)

    # If the change in estimated vector is less than the tolerance level
    # we count it as a violation. We only allow a limited number of such
    # violations.
    if(norm(xhat - xhatprev, 1) < options.tolerance)
      numBelowTolerance += 1
      if numBelowTolerance > options.maxToleranceViolations
        break
      end
    end

  end



end
