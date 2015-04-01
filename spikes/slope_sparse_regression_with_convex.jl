# Solves sparse regression problems with the SLOPE method of Bogdan et al [Bogdan2014].
# Since the prox step in solving SLOPE is a quadratic program we here use Convex.jl
# to solve it. This should in general be an inferior solution to the FastProxSL1 method
# described in [Bogdan2014]. But it is fun and simple and still performs well. :)

# Beware! This is not yet bug-free...

# The key step in both the proximal gradient and the accelerated proximal gradient algorithms
# to solve SLOPE is to solve the prox problem. According to 2.2 and 2.3 in [Bogdan2014] the prox
# can be solved by a quadratic program which we implement in Convex.jl:
using Convex

# Set up x once and for all so we need not do it on each prox optimization round.
setup_x_and_constraints(p) = begin
  x = Variable(p)
  constraints = Constraint[]
  for i in 1:(p-1)
    constraints += (x[i] >= x[i+1])
  end
  constraints += (x[p] >= 0)
  return x, constraints
end

# The prox problem we solve assumes y are sorted and positive. We keep the permutation vector
# around so we can revert later.
sort_and_normalize(y) = begin
  absy = abs(y)
  return sortperm(absy, rev = true), sign(y), absy
end

# Go back from a sortperm to orig order.
function revert_sortperm(y, yperm)
  p = length(y)
  yr = zeros(p)
  for i in 1:p
    yi = yperm[i]
    yr[yi] = y[i]
  end
  yr
end

solve_prox_problem_with_convexjl(yproxvals, lambdas, x, constraints) = begin
  yperm, sngy, absy = sort_and_normalize(yproxvals)
  ynormalized = absy[yperm]

  proxproblem = minimize(0.5 * norm(ynormalized .- x, 2)^2 + sum(lambdas .* x), constraints)
  solve!(proxproblem)
  if proxproblem.status != :Optimal
    throw(string(proxproblem))
  end

  # Reorder and put back the signs...
  sngy .* revert_sortperm(x.value, yperm)
end

# Use a given lambdas sequence or set it up based on the simple sequence from
# Candes2015 paper if not given.
function check_and_setup_lambdasequence(p, lambdas = nothing)
  if lambdas != nothing
    @assert p == length(lambdas)
  else
    # Page 7 of the latest Candes2015 paper on SLOPE states that a very simple lambda sequence can be used:
    lambdafunc(k, p) = sqrt(2*log(p/k))
    lambdas = map(k -> lambdafunc(k, p), 1:p)
  end
  return lambdas
end

# The simplest SLOPE solver is just a proximal gradient. It assumes we have a way to solve the prox
# problem.
function solve_SLOPE_w_proximal_gradient(X::Matrix{Float64}, y::Vector{Float64}, proxsolver::Function;
  stepsizefunc = nothing,       # Default is to use 1/norm(X)^2, so a fixed scheme...
  lambdas = nothing,            # Default is to use the simple sequence below from Candes2015 paper...
  maxiterations = int(5e3), stopval = 1e-3, verbose = true)

  # Assert valid inputs, setup and precalc
  n, p = size(X)
  @assert n == length(y)
  lambdas = check_and_setup_lambdasequence(p, lambdas)
  Xprim = X'        # Precalc for speed 
  bprev = randn(p)  # Initial guess is random
  if stepsizefunc == nothing
    limitval = (2 / norm(X)^2) / 2 # Divide by 2 to ensure it is less than the limitval
    stepsizefunc = (i,b) -> limitval
  end

  # Now iterate for maxiterations steps or until stopval reached.
  for k in 1:maxiterations
    tk = stepsizefunc(k, bprev)
    yproxvals = bprev .- tk * Xprim * (X * bprev .- y)
    bnew = proxsolver(yproxvals, lambdas)
    delta = norm(bnew .- bprev)
    println("k = ", k, ", delta = ", delta)
    if delta < stopval
      return bnew
    end
    bprev = bnew
  end
  return bnew
end

function solve_SLOPE_w_convexjl(X::Matrix{Float64}, y::Vector{Float64}; options...)
  p = size(X, 2)
  x, constraints = setup_x_and_constraints(p)
  proxsolver(yvals, lambdas) = solve_prox_problem_with_convexjl(yvals, lambdas, x, constraints)
  solve_SLOPE_w_proximal_gradient(X, y, proxsolver; options...)
end

# Now generate a random, sparse regression problem
function rand_sparse_problem(p; a = iceil(log(p)), n = int(p/2), amplitude = 10.0, sigma = amplitude/100, shuffle = true)
  beta = vcat(amplitude*ones(a)+randn(a), zeros(p-a))
  indices = collect(1:p)
  if shuffle
    shuffle!(indices)
    beta = beta[indices]
  end
  X = randn(n, p)
  errors = sigma*randn(n)
  y = X * beta .+ errors
  return beta, X, y, indices, n, p, a, sigma, amplitude, errors
end

beta, X, y, indices, n, p, a, sigma, amplitude, errors = rand_sparse_problem(10);

# And let SLOPE loose on it...
@time betahat = solve_SLOPE_w_convexjl(X, y)
println(beta)
println(betahat)
println(norm(beta .- betahat))

# There are bugs though...
