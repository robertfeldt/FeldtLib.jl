# This implements the L0 EM algorithm for regularized regression as described in:
#
#  Zhenqiu Liu and Gang Li. "Efficient Regularized Regression for Variable 
#  Selection with L0 Penalty" Submitted to arXiv on July 29th 2014.
#
# theta = l0EM(X, y, lambda, epsilon, deltaTreshold)
#
# Find coefficients 'theta' for regularized regression with L0 penalty, i.e.
# theta minimizes 0.5*norm(y - X*theta) + lambda/2*sum(theta > 0.0).
#
# Outputs:
#  'theta' is an M*1 vector with the (sparse) coefficients.
#
# Inputs:
#  'X' is an N by M matrix.
#
#  'y' is an N by 1 array.
#
#  'lambda' is the regularization parameter. By default it is set via AIC, i.e. to 2.
#
#  'epsilon' is the minimum coefficient value that is allowed, below which the 
#    coefficient is set to zero. Defaults to 1e-4.
#
#  'delta_treshold' is the treshold value below which the iteration is deemed to have 
#    converged. Defaults to 1e-6.
#
#  'nonnegative' is true iff only non-negative coefficients are allowed. Defaults to false.
function l0EM(X, y, lambda = 2; epsilon = 1e-4, deltaTreshold = 1e-5, nonnegative = false)
  N, M = size(X)
  lambda_eye = lambda * eye(N)
  xt = X'
  theta = xt * inv(X * xt .+ lambda_eye) * y

  if nonnegative
    set_negative_to_zero!(theta)
  end

  while true
    # E-step:
    eta = theta

    # M-step:
    eta_squared = eta .^ 2
    xt_eta = broadcast(*, eta_squared, xt)
    theta = xt_eta * inv(X * xt_eta .+ lambda_eye) * y

    # Ensure non-negativity if required
    if nonnegative
      set_negative_to_zero!(theta)
    end

    # We have converged if change is too small => break.
    if norm(theta - eta) < deltaTreshold
      break
    end
  end

  theta[theta .< epsilon] = 0.0

  theta
end

function set_negative_to_zero!(v)
  v[v .< 0.0] = 0.0
  v
end

# Given mild constraints on the co-linearity of features there is 
# a max lambda value we can calculate.
function find_max_lambda(X, y)
  lambda_from_col_j(j) = (sum(X[:,j] .* y)^2) / (4 * sum(X[:,j].^2))
  N, M = size(X)
  maximum(map(lambda_from_col_j, 1:M))
end

num_selected(v) = sum(v .> 0.0)

# Instead of checking a fixed number of lambda values in a range from 1e-8
# to max_lambda we do a binary search to identify lambda values where there
# is a transition in the number of selected variables.
#
#  cs, ns = adaptive_lambda_regularized_regression(X, y, regressor)
#
# Outputs:
#  'cs' is a dict mapping the tried lambda values to the regressed coefficients for
#    each lambda.
#
#  'ns' is a dict mapping an integer number of selected variables to the lambda value
#    for which that number of selected variables was first encountered.
#
# Inputs:
#  'X' is an N by M matrix corresponding to N cases each with M measured features.
#
#  'y' is an N by 1 array corresponding to the dependent variable for the N cases.
#
#  'regressor' is a Julia function to call to perform regularized regression.
#    It needs to have an interface like regressor(X, y, lambda; kws...).
#    Defaults to using the L0 EM algorithm for regularized regression.
#
# Optional keyword arguments:
#  'minLambda' is the minimum lambda value to use. Defaults to 1e-7 and might not
#    converge at all if a smaller value is selected.
#
#  'numLambdas' is an indication (not exact) of the number of calls to the regressor
#    that we are allowed to do. The actual number of calls can deviate somewhat but not a lot.
#
#  'stepDivisor' determines how small lambda increments are considered when searching in 
#    between the standard log-spaced values of LASSO and similar algorithms. Defaults to
#    2^12 which tends to give a sensible number of alternative models to choose from. 
#    Increase it to find a large number of selected number of variables (returned in 'ns').
#
#  'kws' are any other keyword arguments that will be passed on to the regressor.
#
function adaptive_lambda_regularized_regression(X, y, regressor = l0EM; 
  minLambda = 1e-7, numLambdas = 100, stepDivisor = 2^12,
  kws...)

  max_lambda = find_max_lambda(X, y)
  log_increment = (log(max_lambda + 1) - log(minLambda + 1)) / (numLambdas-1)

  num_vars_to_lambda = Dict{Int64, Float64}()
  coeffs = Dict{Float64, Matrix{Float64}}()

  num_calls_to_regularized_regressor = 0

  update_if_not_there!(l) = begin
    if haskey(coeffs, l)
      cs = coeffs[l]
    else
      coeffs[l] = cs = regressor(X, y, l; kws...)
      num_calls_to_regularized_regressor += 1
    end
    nvars = num_selected(cs)
    if !haskey(num_vars_to_lambda, nvars)
      num_vars_to_lambda[nvars] = l
    end
    return (cs, nvars)
  end

  # Search the lambda values from minLambda to maxLambda
  binary_index_search(minLambda, maxLambda, minDelta) = begin
    #println("bis($minLambda, $maxLambda, $minDelta)")
    if num_calls_to_regularized_regressor <= numLambdas && (maxLambda - minLambda > minDelta)
      csmin, nvarmin = update_if_not_there!(minLambda)
      #println("  nvarmin = $nvarmin ($minLambda)")
      csmax, nvarmax = update_if_not_there!(maxLambda)
      #println("  nvarmax = $nvarmax ($maxLambda)")
      if (nvarmin - nvarmax > 1)
        mid = minLambda + (maxLambda - minLambda) / 2
        csmid, nvarmid = update_if_not_there!(mid)
        #println("  nvarmid = $nvarmid ($mid)")
        if nvarmin - nvarmid > 1
          binary_index_search(minLambda, mid, minDelta)
        end
        if nvarmid - nvarmax > 1
          binary_index_search(mid, maxLambda, minDelta)
        end
      end
    end
  end

  binary_index_search(minLambda, max_lambda, log_increment/stepDivisor)

  return (coeffs, num_vars_to_lambda)
end


N = 50
M = 100
X = randn(N, M)
y = 1.0 * X[:,1] + 2.0 * X[:,2] + 3.0 * X[:,3] + 4.0 * X[:,4] + 0.01 * rand(N, 1)

# 0.25, 0.50, AIC, BIC, RIC, maxLambda
#lambdas = [0.25, 0.50, 2.0, log(N), 2*log(M), find_max_lambda(X,y)]

#ts = map((l) -> l0EM(X, y, l), lambdas)
#t0 = l0EM(X, y, 1e-4)

#function log_split_lambdas(X, y, numLambdas = 10)
#  exp(linspace(0.0, log(find_max_lambda(X, y) + 1), numLambdas)) - 1.0 + 1e-6
#end

#ts100 = map((l) -> l0EM(X, y, l), log_split_lambdas(X, y, 100))
#num_vars = map((t) -> sum(t .> 0.0), ts100)

@time cs, ns = adaptive_lambda_regularized_regression(X, y, l0EM);
#sort(collect(keys(ns)))
#length(cs)
#@time cs, ns8 = adaptive_lambda_regularized_regression(X, y, l0EM; stepDivisor = 2^8);
#@time cs, ns12 = adaptive_lambda_regularized_regression(X, y, l0EM; stepDivisor = 2^12);
#@time cs, ns14 = adaptive_lambda_regularized_regression(X, y, l0EM; stepDivisor = 2^14);
#@time cs, ns16 = adaptive_lambda_regularized_regression(X, y, l0EM; stepDivisor = 2^16);
