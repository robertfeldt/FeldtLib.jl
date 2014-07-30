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
function l0EM(X, y, lambda = 2; epsilon = 1e-4, deltaTreshold = 1e-6, nonnegative = false)
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

# Given mild constraints on the linearity of features there is 
# a max lambda value we can calculate.
function find_max_lambda(X, y)
  lambda_from_col_j(j) = (sum(X[:,j] .* y)^2) / (4 * sum(X[:,j].^2))
  N, M = size(X)
  maximum(map(lambda_from_col_j, 1:M))
end

N = 50
M = 50
X = rand(N, M)
y = 1.0 * X[:,1] + 2.0 * X[:,2] + 5.0 * X[:,3] + 0.01 * rand(N, 1)

# 0.25, 0.50, AIC, BIC, RIC, maxLambda
lambdas = [0.25, 0.50, 2.0, log(N), 2*log(M), find_max_lambda(X,y)]

ts = map((l) -> l0EM(X, y, l), lambdas)

log_split_lambdas(X, y, numLambdas = 10) = exp(linspace(0.0, log(find_max_lambda(X, y) + 1), numLambdas)) - 1.0

# Does not converge!
#ts100 = map((l) -> l0EM(X, y, l), log_split_lambdas(X, y, 10))
