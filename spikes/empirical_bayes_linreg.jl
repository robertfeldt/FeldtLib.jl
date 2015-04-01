# This is a Julia implementation of Empirical Bayes model selection and estimation
# in sparse high-dimensional Gaussian linear regression models as described in the
# paper:
#    Martin, Mess, and Walker, "Empirical Bayes posterior concentration in sparse
#    high-dimensional linear models," arXiv:1406:7718
#

# We need a lasso solution to start from so use GLMNet.
using GLMNet

# Get the lasso solution from a cross-validation run. Since it is only for starting
# the MH chain in empirical bayes below we only run a 5-fold CV.
function lasso_cv(x, y)
  cv = glmnetcv(x, y)
  cv.path.betas[:, indmin(cv.meanloss)]
end

# Prior for complexity/size of a model.
function dcomplex(x, n, p, a, b, calclog = true)
  o = -x * (log(b) + a * log(p)) + log(float(x .<= n))
  return (calclog ? o : exp(o))
end

# Empirical Bayes model selection and estimation Run Metropolis-Hastings chain to sample 
# from marginal posterior.
#
#   S, s, i, beta = empirical_bayes_linreg(x, y; options...)
#
# Inputs:
#  x = predictor variables, n*p, with n samples of p variables
#  y = response variables, n*1, with n samples of response
#
# Outputs:
#
# Options:
#  sigma = error variance, if nothing then it is estimated from data
#  alpha = likelihood fraction (0 < alpha < 1)
#  gamma = conditional prior parameter (0 <= gamma, should be close to 0)
#  numSamples = number of Monte Carlo samples
#  burnIn = fraction of additional burn in samples taken first but not saved
#  sampleCoefs = true iff samples of the coefficients should be returned
#  logF = function to calc log of prior for the model size
#
function empirical_bayes_linreg(X, y; 
  sigma2 = nothing,
  alpha = 0.999,
  gamma = 0.001,
  numSamples = 10000,
  burnIn = 0.25,
  logF = dcomplex
)

  # Check inputs and prep.
  n, p = size(X)
  @assert size(y, 1) == n
  @assert size(y, 2) == 1

  numBurnin = iround(burnIn * numSamples)

  # We will start the chain from the cross-validated Lasso solution so get that.
  lasso_coefs = lasso_cv(X, y)

  if sigma2 == nothing
    z = y - X * lasso_coefs
    sigma2 = sum(z .^ 2) / max(1, n - sum(S))
  end

  twosigma2 = 2 * sigma2
  v = gamma + alpha / sigma2
  logvhalf = log(v) / 2

  # We will use an array of active coefficients to represent a model, S (selected).
  S = find(lasso_coefs .!= 0.0)

  # Get the log priors so we need not recalc them all the time
  logpriors = map((len) -> logF(len) - lchoose(p, len), 0:n)

  # Save the models and which model was used at each step
  models = Any[]        # Each model is unique here and saved as (model, length, betahat, sse, U)
  modelhash = Dict{Vector{Int64}, Int64}() # Map from a model to its index in models array
  sampledmodels = zeros(Int64, numSamples) # Index of model sampled in each step of the chain

  addmodel(model) = begin
    if haskey(modelhash, model)
      entry = models[modelhash[model]]
    else
      # Do linear regression on predictor variables selected by the model
      betas, residuals, sse = linreg_for_subset_of_columns(X, y, model)

      # Create entry for models array
      i = 1 + length(models) # will be at next free index
      modelhash[model] = i
      s = length(model)
      logposterior = logpriors[s] - alpha * sse / twosigma2 - s * logvhalf
      entry = (i, s, model, betas, sse, logposterior)

      # Save entry in models array
      push!(models, entry)
    end
    return entry
  end

  # Add our first model, the one from the Lasso.
  i, s, S, coefs, sse, logposterior = addmodel(S)

  # Now take the Monte Carlo samples
  for j in 1:(numBurnin + numSamples)

    # Get a new model proposal based on random sampling from current model
    Snew = random_proposal_model(S, n, p)

    # Add the model after calculating its implications/variables.
    inew, snew, Snew, coefs, sse, logposteriornew = addmodel(Snew)

    println("model ", j, ": ", snew, ", ", sse, ", ", logposteriornew, "\n")

    # Metropolis-Hastings accept step: Randomly accept proposal dependent on 
    # its posterior compared to the posterior of the current model.
    if rand() <= exp(logposteriornew - logposterior)
      i = inew
      S = Snew
      logposterior = logposteriornew
    end

    # Save sample if we are past the burnin
    if j > numBurnin
      sampledmodels[j - numBurnin] = i
    end

  end

  return models, sampledmodels

end

bigfloat_to_float(bf, digits = 7) = int(bf*(10^digits))/(10^digits)
lchoose(n, m) = bigfloat_to_float(log(binomial(BigInt(n), BigInt(m))))

function linreg_for_subset_of_columns(X, y, model)
  xselected = X[:, model]
  betas = linreg(xselected, y)
  residuals = y .- (betas[1] + xselected * betas[2:end])
  sse = sum( residuals.^2 )
  return betas, residuals, sse
end

# Helper functions for proposal function
delete_at(a, i) = [a[1:(i-1)], a[(i+1):end]]
sample(ary) = ary[rand(1:length(ary))]
add_new_from_range(a, maxval, minval = 1) = push!(copy(a), sample(collect(setdiff(Set(minval:maxval), Set(a)))))

# Proposal function that randomly samples a new model close to given one.
function random_proposal_model(model, maxlen, numvars)
  s = length(model)
  if s == maxlen # If we are at the max length we can only take away variables
    return delete_at(model, rand(1:s))
  elseif s == 1 # If we are at length 1 we can only add
    return sort(add_new_from_range(model, numvars))
  else
    if rand() <= 0.5
      return delete_at(model, rand(1:s))
    else
      return sort(add_new_from_range(model, numvars))
    end
  end
end

# To sample coefficients given a model with a restricted set of
# predictor variables we must first calculate the U as defined 
# in the paper.
function calc_U(xselected, model)
  xtxselected = xselected' * xselected
  V = pinv(xtxselected)
  U = chol(V)
  return U
end

sample_coefficients(s, U, sqrtv, betas) = betas + sqrtv * U * randn(s)
