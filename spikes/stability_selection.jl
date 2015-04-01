using GLMNet

function logrange(min, max, num)
  shift = 1.0 - min
  exp(linrange(log(min + shift), log(max + shift), num)) - shift
end

macro validate_and_init_for_stability_selection()
  esc(quote
    @assert 0.0 <= alpha <= 1.0
    X = float64(X)
    y = float64(y)
    N, P = size(X)
    halfsize = ifloor(N / 2)
    counts = zeros(Int64, P, nlambda)
    updatecounts = zeros(Int64, nlambda)
    Xn, xcolnorms = normalize_columns(X)
    Xnw = zeros(Float64, N, P)
    indices = collect(1:N)
  end)
end

function normalize_columns(A)
  acolnorms = map((j) -> norm(A[:,j]), 1:size(A, 2))
  An = broadcast(/, A, acolnorms')
  return (An, acolnorms)
end

macro update_counts(counts, updatecounts, r)
  quote
    for li in 1:length($(r).lambda)
      $(updatecounts)[li] += 1
      for bi in 1:P
        if $(r).betas[bi, li] != 0.0
          $(counts)[bi, li] += 1
        end
      end
    end
  end
end

# Stability selection in the spirit of Meinshausen&Buhlman but with the Complementary Pairs
# variant described in
#   Rajen Shah and Richard J. Samworth, "Variable selection with error control: Another
#   look at Stability Selection", arxiv.org, 2011, http://arxiv.org/pdf/1105.5578.pdf
#
#   selfreqs = stability_selection(X, y; options...)
#
# Inputs:
#  X is the N*P design matrix
#  y is the N*1 variable to predict
#
# Output:
#  selfreqs is the probability that a feature was selected during the first
#             maxit iterations of the GLMNet LASSO algorithm. The features are
#             randomly reweighted in [alpha, 1.0] and the probability is estimated
#             by numbootstraps samples.
#
# Options:
#  maxit is the maximum number of iterations of the LASSO
#  alpha is the minimum weight of a feature
#  numbootstraps is the number of bootstrap samples
#  parallel is true if it should run the
function stability_selection(X::AbstractMatrix, y::AbstractVector;
                             numbootstraps = 50,
                             maxit = 500,
                             alpha = 0.20,
                             nlambda = 25,
                             lambdas = [],
                             options...)

  @validate_and_init_for_stability_selection

  # Run one glmnet fit with only a few lambdas so we can extract their min and max.
  if length(lambdas) == 0
    res = glmnet(Xn, y; nlambda = 3, options...)
    lambdamax, lambdamin = res.lambda[1], res.lambda[end]
    lambdas = reverse(logrange(lambdamin, lambdamax, nlambda))
  end

  # Now take numbootstraps (parallel) samples
  map(1:numbootstraps) do i

    # Reweight each covariate in each round. We can do this via the penalty_factors
    # option instead?!
    weights = alpha + (1.0 - alpha) * rand(P)
    for j in 1:P
      for i in 1:N
        Xnw[i, j] = Xn[i, j] / weights[j]
      end
    end

    # Randomly split sample in two sets
    shuffle!(indices)
    idxs1 = indices[1:halfsize]
    idxs2 = indices[(halfsize+1):(2*halfsize)] # Can't take to end since they need to be same size

    # Run GLMNet on each of the samples and count which variables are selected
    r1 = glmnet(Xnw[idxs1,:], y[idxs1]; lambda = lambdas, standardize = false, intercept = true, maxit = maxit)
    r2 = glmnet(Xnw[idxs2,:], y[idxs2]; lambda = lambdas, standardize = false, intercept = true, maxit = maxit)

    @update_counts(counts, updatecounts, r1)
    @update_counts(counts, updatecounts, r2)
  end

  # Complementary Pairs SS just sums all counts per feature
  #freqs = sum(counts, 2) / sum(updatecounts)

  # Original muhleinbein paper used the max??
  freqs = broadcast(*, counts, (1 ./ updatecounts)')
  # Return the max frequency achieved per feature
  freqs = maximum(freqs, 2)

  freqs[:]
end

# This implements the treshold (tau) selection procedure as described
# at the end of page 10 in the Rajan&Shamworth paper.
#
# Input:
#  P is the number of original features
#  l is max number of low probability features selected, on average
#  B is number of bootstraps, i.e. 2*B is the number of subsets used
function select_treshold(P, l = 0.50, B = 50)
  q = sqrt(0.8*l*P)
  minD(q/P, B)
end

function minD(theta, B, rs = [-0.5, -0.25])
  minimum(vcat(ones(B), tailprobabilities(theta^2, B, rs[1]), tailprobabilities(theta, B, rs[2])))
end

# Calculate the tail probabilities (based on r-concave funcs) for each
# relevant tau. Based on Shah and Shamworth paper and their R code:
#   http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
function tailprobabilities(eta, B, r)
  MAXa = 100000
  MINa = 0.00001
  s = -1/r
  etaB = eta * B
  k_start = iceil(2 * etaB) + 1
  output = ones(B)
  if k_start > B
    return(output)
  end

  a_vec = MAXa * ones(B)


end

function indover(v, treshold = 0.70)
  collect(1:length(v))[v .> treshold]
end

N, P, k = (50, 200, 5)
coefs = shuffle(vcat(100*rand(k), zeros(P-k)))
X = randn(N, P)
errors = 1*randn(N)
y = X * coefs .+ errors

tic()
selfreqs = stability_selection(X, y; maxit = 500, alpha = 0.20, nlambda = 20)
t = toq()
reverse(sort(coefs))[1:k]
reverse(sort(selfreqs))[1:(2k)]
indover(coefs, 0.0)
indover(selfreqs, min(N/P, 1.0)) # Shah&Samworth recommends N/P as the treshold, I modified to min(N/P, 1.0) since frequencies <= 1.0
