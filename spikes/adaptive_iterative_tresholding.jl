# This implements the Adaptive Iterative Tresholding algorithm for
# compressed sensing / regularized regression. This implementation
# is based on the paper:
#  Yu Wang, Jinshan Zeng, Zhimin Peng, Xiangyu Chang, Zongben Xu,
#  "On Linear Convergence of Adaptively Iterative Thresholding 
#  Algorithms for Compressed Sensing", arxiv, 2014, http://arxiv.org/abs/1408.6890
#
function adaptive_iterative_tresholding(A, b; kws...)
  M, N, k, xprev, x, An, acolnorms, sAnt, f, maxiterations, tolerance = unpack_ait_parameters(A, b; kws...)

  # Since we need one iteration before we can judge convergence
  # we check for convergence last in the loop.
  t = 0
  while true
    t += 1
    z = xprev .+ sAnt * (b .- An * xprev)
    zp = sortperm(abs(z), rev = true)
    tau = abs(z[zp[k+1]])

    # Treshold x
    for i in 1:N
      ip = zp[i]
      if i <= k
        x[ip] = f(tau, z[ip])
      else
        x[ip] = 0.0
      end
    end

    # Check convergence
    diff = norm(x .- xprev)
    if t >= maxiterations || diff < tolerance
      return (x ./ acolnorms, t)
    end

    # Switch their meaning => we need not allocate new memory in each iteration
    xprev, x = x, xprev
  end
end

function unpack_ait_parameters(A, b;
  k = -1,             # Sparsity level, this should be set as an upper bound on the number of non-zero coefficients of x. If set below the actual sparseness level the algorithm might not find right solution.
  s = 1.0,            # Step length, should be ~ in range [0.4, 1.3] but 1.0 seem to be generally good
  x0 = nothing,       # A starting point, will be randomized if not given
  f::Function = scad, # Defining function used in tresholding, by default SCAD
  tolerance = 1e-6,   # Convergence level, i.e. terminate algorithm if norm change is less than this between iterations
  maxIterations = 1e3 # Maximum number of iterations allowed. Often it converges in ~20-30 iterations if ever.
  )

  # Check and then initialize all parameters. Note that we use N for the num of covariates and
  # M for the number of cases. This is in line with the paper above but it is more common to use
  # N = M and P = N.
  M, N = size(A)
  @assert length(b) == M
  k = (k == -1) ? ifloor(sqrt(N)) : k
  @assert 1 <= k < N
  xprev = (x0 == nothing) ? zeros(N) : copy(x0)
  x = zeros(Float64, N)

  An, acolnorms = normalize_columns(A)
  sAnt = s * An'

  return (M, N, k, xprev, x, An, acolnorms, sAnt, f, maxIterations, tolerance)
end

function normalize_columns(A)
  acolnorms = map((j) -> norm(A[:,j]), 1:size(A, 2))
  An = broadcast(/, A, acolnorms')
  return (An, acolnorms)
end

# Inner AIT iteration loop. Assumes An has been normalized.
function ait_iteration_loop!(An, b, k, N, M, sAnt, tolerance, maxiterations, xprev, trainindices, Ntrain, f)
  t = 0
  x = zeros(Float64, N)

  # For now, later we will loop over the trainindices
  btrain = b[trainindices]
  Antrain = An[trainindices, :]
  sAnttrain = sAnt[:, trainindices]

  while true
    t += 1
    z = xprev .+ sAnttrain * (btrain .- Antrain * xprev)
    zp = sortperm(abs(z), rev = true)
    tau = abs(z[zp[k+1]])

    # Treshold x
    for i in 1:N
      ip = zp[i]
      if i <= k
        x[ip]::Float64 = f(tau::Float64, z[ip]::Float64)
      else
        x[ip]::Float64 = 0.0
      end
    end

    diff = norm(x .- xprev)
    if t >= maxiterations || diff < tolerance
      return (x, t)
    end

    xprev, x = x, xprev
  end
end

function ait2(A, b; kws...)
  M, N, k, xprev, x, An, acolnorms, sAnt, f, maxiterations, tolerance = unpack_ait_parameters(A, b; kws...)
  x, its = ait_iteration_loop!(An, b, k, N, M, sAnt, tolerance, maxiterations, xprev, collect(1:M), M, f)
  return (x ./ acolnorms, its)
end

function scad_penalty_tresholding(a, tau, u)
  absu = abs(u)
  if absu > a * tau
    u
  elseif absu > 2 * tau
    ((a-1)*u - sign(u)*a*tau)/(a-2)
  elseif absu > tau
    absu - sign(u)*tau
  else
    0.0
  end
end

scad(tau, u) = scad_penalty_tresholding(3.7, tau, u)
hard(tau, u) = abs(u) > tau ? u : 0.0

splitindices(n, ratio = 0.10) = begin
  randindices = shuffle(collect(1:n))
  numtest = ifloor(ratio*n)
  test = randindices[1:numtest]
  train = randindices[(numtest+1):end]
  return (train, test)
end

# Run AIt for multiple k-values from 2*sqrt(P) down to 1
# and then select the one with lowest MSE on a TestRatio
# holdout test set.
function ait_best_k(A, b, TestRatio = 0.20; kws...)
  N, P = size(A)

  # Split in train and test sets
  trainidxs, testidxs = splitindices(N, TestRatio)
  Atrain = A[trainidxs,:]
  btrain = b[trainidxs]
  Atest = A[testidxs,:]
  btest = b[testidxs]

  maxP = ifloor(2*sqrt(P))
  if maxP >= P
    maxP = P - 1
  end

  local bestcoefs = bestits = bestmse = Inf

  mses = Dict{Int64,Float64}()

  for k in maxP:-1:1
    coefs, its = adaptive_iterative_tresholding(Atrain, btrain; k = k, kws...)
    mses[k] = mse(btest, Atest * coefs)
    if mses[k] < bestmse
      bestcoefs = coefs
      bestits = its
      bestmse = mses[k]
    end
  end

  return bestcoefs, bestits, mses
end

# using PyPlot
# x = linspace(-5.0, 5.0, 200)
# y = map((x) -> scad_penalty_tresholding(3.7, 1.0, x), x)
# plot(x, y, color="red", linewidth=2.0, linestyle="--")

mse(y,yhat) = mean( (y .- yhat).^2 )

# N-fold cross-validation to select the sparseness level for AIT.
function nfoldCV_ait(A, b)
  M, N = size(A)

  num_ks = 5
  ks = ifloor(exp(log(N-1) / num_ks * reverse(1:num_ks)))
  xs = zeros(num_ks, N)
  x0 = zeros(N)
  its = zeros(Int64, num_ks)

  for i in 1:num_ks
    x0, its[i] = adaptive_iterative_tresholding(A, b; k=k, x0 = x0)
    xs[i,:] = x
  end
end

function ait_mse(A, b, k)
  betas, its = adaptive_iterative_tresholding(A, b; k = k)
  return (betas, mse(b, A * betas))
end

function adaptive_k_ait(A, b; kws...)
  N, P = size(A)

  mses = Dict{Int64, Float64}()
  betas = Dict{Int64, Vector{Float64}}()
  mse_for_k(k) = begin
    if !haskey(mses, k)
      betas[k], mses[k] = ait_mse(A, b, k)
    end
    mses[k]
  end 

  binsearch(lowk, highk) = begin
    mselow = mse_for_k(lowk)
    msehigh = mse_for_k(highk)
    if mselow > 1e3 * msehigh
      if highk == lowk + 1
        return (highk, msehigh)
      end
      midk = lowk + ifloor((highk-lowk)/2)
      k, mseb = binsearch(lowk, midk)
      if k == -1
        return binsearch(midk, highk)
      else
        return (k, mseb)
      end
    else
      return (-1, nothing)
    end
  end

  k, mse = binsearch(1, ifloor(min(P/2, 2*sqrt(P))))

  return (k, mse, betas[k], mses)
end

function test_ait_on_random_problem(N, P, kstar)
  A = randn(N, P)
  nzw = map((i) -> rand(10.0:100.0) * (rand() < 0.5 ? 1.0 : (-1.0)), 1:kstar) 
  betas = shuffle([nzw, zeros(P-kstar)])
  b = A * betas

  k, mse, betasp, mses = adaptive_k_ait(A, b)

  show(mses)
  println()

  return (N, P, kstar, k, mse, coeffquality(betas, betasp))
end

#for kstar in [3, 4, 5, 6]
#  for P in [10, 20, 50, 100]
#    println(test_ait_on_random_problem(100, P, kstar), "\n")
#  end
#end

# N, P = (rand(20:100), 20)
# kstar = ifloor(sqrt(P)) - 1
# A = randn(N, P)
# nzw = rand(1.0:10.0) * rand(-1.0:0.01:1.0, kstar) 
# betas = shuffle([nzw, zeros(P-kstar)])
# b = A * betas
# k = -1.0
# s = 1.0
# x0 = nothing
# f = hard
# tolerance = 1e-5
# maxIterations = 1e4
# 
# # Ensure compiled
# adaptive_iterative_tresholding(A[1:10,:], b[1:10])
# ait2(A[1:10,:], b[1:10])
# 
# tic()
# betasp, its = adaptive_iterative_tresholding(A, b)
# t1 = toq()
# tic()
# betasp2, its2 = ait2(A, b)
# t2 = toq()
# nzwp = betasp[abs(betasp) .> 0.0]
# bhat = X' * betasp
# bhat2 = X' * betasp2
# 
# 
# show( [betas betasp betasp2] )
# println("\nN = ", N, ", P = ", P)
# println("MSE1: ", mse(b, bhat))
# println("MSE2: ", mse(b, bhat2))
# println("b1: ", coeffquality(betas, betasp))
# println("b2: ", coeffquality(betas, betasp2))
# println("Timediff: ", t1-t2)