# Info object return after a SLOPE solve.
type SLOPESolveInfo
  runtime       # Runtime
  Aprods        # Number of products with A
  ATprods       # Number of products with A^T
  objPrimal     # Primal objective
  objDual  # Dual objective (possibly for infeasible dual point)
  infeas        # Dual infeasibility
  status        # Status: 1 = optimal, 2 = iterations
end

# Constants for exit status
const STATUS_RUNNING   = 0
const STATUS_OPTIMAL   = 1
const STATUS_ITERATIONS = 2
const STATUS_MSG = ["Optimal", "Iteration limit reached"]

# The default choice for lambdas is Benjamini-Hochberg sequence where q is the FDR:
using Distributions
NormalDistr = Normal()
bh_lambdafunc(i, p, q = 0.10) = quantile(NormalDistr, 1 - i*q/(2*p))

# Lets implement the stack-based FastProx SL1 alg of Bogdan et al to see if that helps...
# This is Alg 4 in the Bogdan2014 paper.
function fast_prox_sl1(yproxvals, lambdas)
  p = length(yproxvals)
  tuples = Any[]
  t = 0

  # Find optimal group levels
  for k in 1:p
    t += 1
    i, j, s = k, k, (yproxvals[k] - lambdas[k])
    w = abs(s)
    push!(tuples, (i, j, s, w))
    while (t > 1) && (tuples[t-1][4] .<= w)
      itm1, jtm1, stm1, wtm1 = tuples[t-1]
      it, jt, st, wt = tuples[t]
      wnew = (jtm1 - itm1 + 1)/(jt - itm1 + 1) * stm1 .+ (jt - it + 1)/(jt - itm1 + 1) * st
      tuples[t-1] = (itm1, jt, stm1 .+ st, abs(wnew))
      pop!(tuples) # Delete last entry
      t -= 1
    end
  end

  # Set entries in x for each block
  x = zeros(p)
  for l in 1:t
    for k in tuples[l][1]:tuples[l][2]
      x[k] = tuples[l][4]
    end
  end
  x
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

solve_prox_problem_with_fastproxsl1(yproxvals, lambdas) = begin
  yperm, sngy, absy = sort_and_normalize(yproxvals)
  ynormalized = absy[yperm]

  xvals = fast_prox_sl1(ynormalized, lambdas)

  # Reorder and put back the signs...
  sngy .* revert_sortperm(xvals, yperm)
end

# Test values for params
#iterations = 10000      # Maximum number of iterations
#q = 0.05
#verbosity = 1         # 0 = nothing, 1 = major, 2 = every
#optimIter = 1         # Iterations between optimality-condition checks
#gradIter = 20         # Iterations between full gradient computations
#tolInfeas = 1e-6       # Maximum allowed dual infeasibility
#tolRelGap = 1e-6       # Stopping criterion for relative primal-dual gap
#xInit = zeros(size(A,2))  # Initial value of x
#lambdas = nothing

# solve_SLOPE  Sorted L1 parameter estimation solver
#
# x, info = solve_SLOPE(A, b, lambda; options) solves the Slope problem
#
#     Minimize 1/2*||Ax-b||_2^2 + sum_i (lambda_i * |x|_[i])
#
# where |x|_[i] denotes the i-th largest entry in |x|. The entries in
# lambda must be nonnegative and in non-increasing order. When lambda is a
# scalar, the above formulation is reduced to the Lasso:
#
#     Minimize 1/2*||Ax-b||_2^2 + lambda * ||x||_1.
#
# Original MatLab version is:
# Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
# Julia version by Robert Feldt, robert.feldt@gmail.com
# was based on Adlas.m in SLOPE Toolbox version 1.0 while the FastProxSL1 algorithm
# was implemented from scratch in Julia by Robert based on the Bogdan et al 2014 paper 
# on SLOPE. We also excluded the Lasso from the implementation since there are already
# LASSO implementations for Julia.
#
#   The SLOPE Toolbox is free software: you can redistribute it
#   and/or  modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation, either version 3 of
#   the License, or (at your option) any later version.
#
#   The SLOPE Toolbox is distributed in the hope that it will
#   be useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with the SLOPE Toolbox. If not, see
#   <http://www.gnu.org/licenses/>.
#
function solve_SLOPE(A, b, lambdas = nothing; 
  iterations = 10000,      # Maximum number of iterations
  q = 0.05,                 # Expected FDR, False Discovery Rate
  verbosity = 1,         # 0 = nothing, 1 = major, 2 = every
  optimIter = 1,         # Iterations between optimality-condition checks
  gradIter = 20,         # Iterations between full gradient computations
  tolInfeas = 1e-6,       # Maximum allowed dual infeasibility
  tolRelGap = 1e-6,       # Stopping criterion for relative primal-dual gap
  xInit = zeros(size(A,2))  # Initial value of x
  )

# Start timer
t0 = tic()

# Get problem dimension
n = size(A, 2)

# Create lambda sequence unless already done.
if lambdas == nothing
  lambdas = map(k -> bh_lambdafunc(k, n, q), 1:n)
end

# Ensure that lambda is non-increasing, nonnegative and has at least one non-zero entry.
if ((length(lambdas) > 1) && any(lambdas[2:end] .> lambdas[1:end-1]))
  throw("Lambdas must be non-increasing.")
end
if (lambdas[end] < 0)
  throw("Lambdas must be nonnegative")
elseif (lambdas[1] == 0)
  throw("Lambdas must have at least one nonnegative entry.")
end

# Get initial lower bound on the Lipschitz constant
x = randn(n,1)
x = x / norm(x,2)
x = A'*(A*x)
L = norm(x,2)

# Initialize parameters and iterates
t     = 1
eta    = 2
lambdas  = lambdas[:]
b     = b[:]
x     = xInit
y     = x
Ax    = A*x
fPrev  = Inf
iter   = 0

info = SLOPESolveInfo(
  0.0, # runtime is set after finishing
  2,   # Aprods
  1,   # ATprods
  0.0, # objPrimal is set later
  0.0, # info.objDual   is set later
  0.0, # Dual infeasibility
  STATUS_RUNNING
)

solve_prox = solve_prox_problem_with_fastproxsl1

if (verbosity > 0)
  @printf("%5s  %9s  %9s  %9s  %9s", "Iter", "||r||_2", "Gap", "Infeas.", "Rel. gap")
end

# Main loop
while (true)

  # Compute the gradient at f(y)
  if (mod(iter, gradIter) == 0) # Includes first iterations
   r = A*y - b
  else
   r = (Ax + ((tPrev - 1) / t) * (Ax - AxPrev)) - b
  end
  g = A' * r
  f = dot(r, r) / 2
  
  # Increment iteration count
  iter += 1

  # Check optimality conditions
  if ((mod(iter, optimIter) == 0))
    # Compute 'dual', check infeasibility and gap
    gs    = sort(abs(g), rev = true)
    ys    = sort(abs(y), rev = true)
    info.infeas = max(maximum(cumsum(gs-lambdas)),0)
    
    # Compute primal and dual objective
    info.objPrimal = ( f + lambdas'*ys)[1]
    info.objDual   = (-f - r'*b       )[1]
    
    if (verbosity > 0)
      str = @sprintf("  %9.2e  %9.2e  %9.2e", info.objPrimal - info.objDual, info.infeas/lambdas[1], 
        abs(info.objPrimal - info.objDual) / max(1,info.objPrimal))
    end
    
    # Check primal-dual gap
    if ((abs(info.objPrimal - info.objDual)/max(1,info.objPrimal) < tolRelGap) && 
        (info.infeas < tolInfeas * lambdas[1]))
      info.status = STATUS_OPTIMAL;
    end
  else
    str = ""
  end

  if ((verbosity == 2) || ((verbosity == 1) && (mod(iter,optimIter) == 0)))
    @printf("%5d  %9.2e%s\n", iter, f, str)
  end

  # Stopping criteria
  if (info.status == STATUS_RUNNING)
    if (iter >= iterations)
      info.status = STATUS_ITERATIONS
    end
  end
  
  if (info.status != 0)
    if (verbosity > 0)
      @printf("Exiting with status %d -- %s\n", info.status, STATUS_MSG[info.status])
    end
    break;
  end
  
  # Keep copies of previous values
  AxPrev = Ax
  xPrev  = x
  fPrev  = f
  tPrev  = t
  
  # Lipschitz search
  while (true)
    # Compute prox mapping
    x = solve_prox(y - (1/L)*g, lambdas/L)
    d = x - y
    
    Ax = A*x
    r  = Ax-b
    f  = dot(r, r)/2
    q  = fPrev + dot(d, g) + (L/2)*dot(d, d)
    
    info.Aprods += 1
    
    if (q >= f*(1-1e-12))
      break
    else
      L = L * eta
    end
  end
  
  # Update
  t = (1 + sqrt(1 + 4*t^2)) / 2
  y = x + ((tPrev - 1) / t) * (x - xPrev)
end

# Final touchup of info values
info.runtime = toc()
info.Aprods = info.Aprods + ceil(iter / gradIter)
info.ATprods  = info.ATprods + iter

# ...and we are finished!
return y, info

end # function solve_SLOPE

function rand_sparse_problem(p; numactive = iceil(log(p)), 
  n = int(p/2), amplitude = 4.1, sigma = 1.0, shuffle = true)

  beta = vcat(amplitude*ones(numactive)+randn(numactive), zeros(p-numactive))
  indices = collect(1:p)
  if shuffle
    shuffle!(indices)
    beta = beta[indices]
  end
  X = randn(n, p)
  errors = sigma*randn(n)
  y = X * beta .+ errors
  return beta, X, y, indices, n, p, numactive, sigma, amplitude, errors
end

# Find the sharpest difference between the found values and use that as the estimate
# of the num active values.
function estimate_active(xhat)
  relsize = abs(xhat) / maximum(abs(xhat))
  relsizeperm = sortperm(relsize, rev=true) 
  srelsize = relsize[relsizeperm]
  relfactors = srelsize[1:end-1] ./ srelsize[2:end]
  numactive = find(rf -> rf == maximum(relfactors), relfactors)[1]
  return relsizeperm[1:numactive]
end

# First estimate the number of active coefficients by a SLOPE solve then
# do a linreg on only those factors.
function sparse_regression_SLOPE(A, b; iterations = 250, options...)
  xhat, info = solve_SLOPE(A, b; iterations = iterations, options...)
  active = estimate_active(xhat)
  Asub = A[:, active]
  res = linreg(Asub, b)
  xhatsub = zeros(size(A, 2))
  xhatsub[active] = res[2:end]
  return xhatsub, res[2:end], active
end

x, A, b, indices, n, p, numactive, sigma, amplitude, errors = rand_sparse_problem(1000; 
  shuffle = false, numactive = 10);

# Split in train and test set
idxs = shuffle(collect(1:n))
trainsize = ifloor(0.80*n)
trainidxs = idxs[1:trainsize]
testidxs = idxs[(trainsize+1):end]
Atrain = A[trainidxs, :]
btrain = b[trainidxs]
Atest = A[testidxs, :]
btest = b[testidxs]

xhat, xselected, selected = sparse_regression_SLOPE(Atrain, btrain)

println(x[1:numactive])
println(xhat[1:numactive])
println("Norm(x .- xhat) = ", norm(x .- xhat))
println("Norm(b .- bhat) = ", norm(btrain .- bhat))
mape(y, yhat) = mean(100.0 * abs(yhat - y)/y)
@printf("Num active = %d, Num active selected = %d", numactive, length(selected))
@printf("Train MAPE = %.2f%%", mape(btrain, Atrain*xhat))
@printf("Test MAPE = %.2f%%", mape(btest, Atest*xhat))
