function forward_stagewise_regression(X, y; epsilon = 0.01, numIterations = 100)
  n, p = size(X)
  ny = length(y)
  @assert n == ny
  forward_stagewise_regression_inner(X, y, n, p, epsilon, numIterations)
end

function forward_stagewise_regression_inner(X, y, n, p, epsilon, numIterations)
  betas = zeros(p, numIterations)
  for k in 2:numIterations
    maxabsinnerproduct = 0.0
    bestinnerproduct = 0.0
    residuals = y .- X * betas[:, k-1]
    bestcolumn = rand(1:p)
    for j in 1:p
      innerproduct = sum( X[:,j] .* residuals )
      if abs(innerproduct) > maxabsinnerproduct
        bestinnerproduct = innerproduct
        maxabsinnerproduct = abs(innerproduct)
        bestcolumn = j
      end
    end
    betas[:, k] = betas[:, k-1]
    betas[bestcolumn, k] += epsilon * sign(bestinnerproduct)
  end
  betas
end

N = 20
P = 5
X = randn(N, P)
y = 1.0 * X[:,1] + 2.0 * X[:,3]
betas = forward_stagewise_regression(X, y; epsilon = 0.10)