# Implements DDS as described in the paper
#  Tolson and Shoemaker, "Dynamically Dimensioned Search Algorithm...", 2007
function dynamically_dimensioned_search(func, minBounds, maxBounds; 
  m = 1000, 
  r = 0.6, # Our experiments shows that often good to start higher with r when gradually changing it
  max_r_factor = 1.0, # Set this higher and the next one lower to scale also the step sizes dynamically
  min_r_factor = 0.01,
  xinitial = nothing,
  sense = :minimize,
  trace = true
  )

  D = length(minBounds)
  @assert D == length(maxBounds)

  deltas = maxBounds .- minBounds

  if isa(xinitial, Nothing)
    x = minBounds .+ rand(D) .* deltas
  else
    @assert length(xinitial) == D
    x = xinitial
  end

  fbest = func(x)

  if sense == :minimize
    fitness_is_better = <
  else
    fitness_is_better = <
  end

  included_vars = ones(D)

  delta_r_factor = max_r_factor - min_r_factor

  for i in 1:m
    # Sample J decision vars to be perturbed
    p = 1 - log(i)/log(m)
    noneincluded = true
    for j in 1:D
      if rand() < p
        included_vars[j] = 1.0
        noneincluded = false
      else
        included_vars[j] = 0.0
      end
    end

    # Ensure at least one selected
    if noneincluded
      included_vars[rand(1:D)] = 1.0
    end

    # Randomly sample and scale to get perturbarion vector
    perturbations = randn(D) .* included_vars .* ((min_r_factor + (p * delta_r_factor)) * deltas)

    # Create new solution vector
    xnew = x .+ perturbations

    # Ensure it is within bounds
    for j in 1:D
      # Ensure it is above min bound
      if xnew[j] < minBounds[j]
        xnew[j] = minBounds[j] + (minBounds[j] - xnew[j])
        if xnew[j] > maxBounds[j]
          xnew[j] = minBounds[j]
        end
      end

      # Ensure it is below max bound
      if xnew[j] > maxBounds[j]
        xnew[j] = maxBounds[j] - (xnew[j] - maxBounds[j])
        if xnew[j] < minBounds[j]
          xnew[j] = maxBounds[j]
        end
      end
    end

    # Evaluate new solution candidate
    fnew = func(xnew)

    # Update if better than previous
    if fitness_is_better(fnew, fbest)
      if trace
        print("$i: found better! fitness = $(round(fnew, 4))")
        if D < 6
          print(" ", xnew)
        end
        println()
      end
      x = xnew
      fbest = fnew
    end
  end

  return (x, fbest)
end

function multi_dds(numReps, func, minBounds, maxBounds; options...)
  fitnesses = zeros(numReps)
  for rep in 1:numReps
    x, fitnesses[rep] = dynamically_dimensioned_search(func, minBounds, maxBounds; options...)
  end
  return mean(fitnesses), std(fitnesses), fitnesses
end

function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end

function ackley(x)
  D = length(x)
  try
    20 - 20.*exp(-0.2.*sqrt(sum(x.^2)/D)) - exp(sum(cos(2 * π * x))/D) + e
  catch e
    # Sometimes we have gotten a DomainError from the cos function so we protect this call
    println(e)
    println("For input x = ", x)
    # Return a very large fitness value to indicate that this is NOT the solution we want.
    # TODO: Fix better way to handle this!
    9.99e100
  end
end

function rosenbrock(x)
  n = length(x)
  return( sum( 100*( x[2:n] - x[1:(n-1)].^2 ).^2 + ( x[1:(n-1)] - 1 ).^2 ) )
end

D = 30
maxBs = 10.0 * ones(D)
minBs = -maxBs
#NumReps = 10
#for N in [5e2, 1e3, 1e4, 1e5]
#  m1, s1, fs1 = multi_dds(NumReps, rastrigin, minBs, maxBs; 
#    m = N, max_r_factor = 1.0, min_r_factor = 1.0, trace = false)
#  m2, s2, fs2 = multi_dds(NumReps, rastrigin, minBs, maxBs; 
#    m = N, max_r_factor = 1.0, min_r_factor = 0.05, trace = false)
#  m3, s3, fs3 = multi_dds(NumReps, rastrigin, minBs, maxBs; 
#    m = N, max_r_factor = 5.0, min_r_factor = 0.05, trace = false)
#
#  println("$(int(N)): default = $(m1) ($(s1)), 1.0-0.05 = $(m2) ($(s2)), 5.0-0.05 = $(m3) ($(s3))")
#end

# Lets use DDS to optimize itself!
function sod(func, D, N, m = 200, NumReps = 10, mode = :minstd, rerunbase = false)
  maxBs = 10.0 * ones(D)
  minBs = -maxBs

  function mdds(x)
    r, maxrf, minrf = x
    m, s, fs = multi_dds(NumReps, func, minBs, maxBs; 
      m = N, r = r, max_r_factor = maxrf, min_r_factor = minrf, trace = false)
    print(".")
    if mode == :minstd
      minf = minimum(fs)
      minf + s
    elseif mode == :meanstd
      m + s
    else
      m
    end
  end

  res = dynamically_dimensioned_search(mdds, [1e-100, 1.0, 1e-100], [1.0, 10.0, 1.0]; m = m, trace = true)

  if rerunbase
    r, maxrf, minrf = res[1]
    println("Rerun with best params:")
    xbest, fbest = dynamically_dimensioned_search(func, minBs, maxBs;
      m = N, r = r, max_r_factor = maxrf, min_r_factor = minrf, trace = true)
    println("Rerun best fitness: $fbest")
  end

  res
end

D = 10
maxBs = 10.0 * ones(D)
minBs = -maxBs
dynamically_dimensioned_search(rosenbrock, minBs, maxBs;
  m = 3e4, r = 0.96, max_r_factor = 6.4, min_r_factor = 3.4, trace = true)

sod(rosenbrock, 10, 1e4, 10, 5, :minstd, true)

# Rastrigin as base problem:
#
# m = 500:
#
# D=  2, N=1e3: ([0.636593,1.06659,0.001],0.0009316716009318879)
# D=  5, N=1e3: ([0.175887,1.09842,0.0194187],0.2693790190666025)
# D= 10, N=1e3: ([0.154905,1.59847,0.001],2.704610222970387)
# D= 20, N=1e3: ([0.22993,1.19872,0.0402445],21.120532067530334)
# D= 30, N=1e3: ([0.335524,1.06239,0.0357023],64.28515634368684)
# D=100, N=1e3: ([0.823432,1.02139,0.0593126],974.3730052913571)
#
# m = 200: 
# D=  2, N=1e4: ([0.855643,1.0,0.00892433],2.4715624989736813e-5)
# D=  5, N=1e3: ([0.075286,1.0395,0.00171595],0.0002987684309726133)
# D= 10, N=1e3: ([0.643196,1.04729,0.0027357],0.002647770416951525)
# D= 20, N=1e3: ([0.972712,1.00394,0.00350361],0.037975116196984224)
# D= 30, N=1e3: ([0.291904,1.28662,0.00880971],1.4776366012816555)
# D=100, N=1e3: ([0.637511,1.18537,0.0610876],127.02431186954585)

# Ackley as base problem:
#
# m = 500:
# D=  2, N=1e3: ([0.609718,1.02179,0.00304914],0.008687833131616746)
# D=  5, N=1e3: ([0.385976,1.07308,0.001],0.019204566095652638)
# D= 10, N=1e3: ([0.612191,1.05159,0.001],0.0668592553112163)
# D= 20, N=1e3: ([0.993921,1.41519,0.001],0.48707509475649696)
# D= 30, N=1e3: ([0.280199,1.20044,0.0131395],1.5370911762677564)
# D=100, N=1e3: ([0.0803656,1.04094,0.0449066],6.985731867803416)
#
# m = 200:
# D=  2, N=1e4: ([0.1506,1.10013,0.00269322],0.0006238838520917218)
# D=  5, N=1e4: ([0.478943,8.7339,0.0096099],0.013788193323540027)
# D= 10, N=1e4: ([0.482863,1.3944,0.00225037],0.005940879821021374)
# D= 20, N=1e4: ([0.0957213,1.07752,0.001],0.008147361361135844)
# D= 30, N=1e4: ([0.538295,1.01522,0.001],0.014400153397059436)
# D=100, N=1e4: ([0.384428,1.03398,0.001],0.12857952728246763)
