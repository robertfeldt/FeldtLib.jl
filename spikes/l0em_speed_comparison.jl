function l0EM_orig(X, y, lambda = 2; 
  epsilon = 1e-4, deltaTreshold = 1e-5, 
  nonnegative = false, maxIterations = 10000)

  N, M = size(X)
  lambda_eye = lambda * eye(N)
  xt = X'
  theta = xt * inv(X * xt .+ lambda_eye) * y

  if nonnegative
    set_negative_to_zero!(theta)
  end

  iterations = 0

  while iterations < maxIterations
    iterations += 1

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

  theta[abs(theta) .< epsilon] = 0.0

  theta
end

function l0EM_us(X, y, lambda = 2; 
  epsilon = 1e-4, deltaTreshold = 1e-5, 
  nonnegative = false, maxIterations = 10000)

  N, M = size(X)
  lambda_eye = UniformScaling(lambda)
  xt = X'
  theta = xt * inv(X * xt .+ lambda_eye) * y

  if nonnegative
    set_negative_to_zero!(theta)
  end

  iterations = 0

  while iterations < maxIterations
    iterations += 1

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

  theta[abs(theta) .< epsilon] = 0.0

  theta
end

function l0EM_solve(X, y, lambda = 2; 
  epsilon = 1e-4, deltaTreshold = 1e-5, 
  nonnegative = false, maxIterations = 10000)

  N, M = size(X)
  lambda_eye = lambda * eye(N)
  xt = X'
  theta = xt * ((X * xt .+ lambda_eye) \ y)

  if nonnegative
    set_negative_to_zero!(theta)
  end

  iterations = 0

  while iterations < maxIterations
    iterations += 1

    # E-step:
    eta = theta

    # M-step:
    eta_squared = eta .^ 2
    xt_eta = broadcast(*, eta_squared, xt)
    theta = xt_eta * ((X * xt_eta .+ lambda_eye) \ y)

    # Ensure non-negativity if required
    if nonnegative
      set_negative_to_zero!(theta)
    end

    # We have converged if change is too small => break.
    if norm(theta - eta) < deltaTreshold
      break
    end
  end

  theta[abs(theta) .< epsilon] = 0.0

  theta
end

function l0EM_us_solve(X, y, lambda = 2; 
  epsilon = 1e-4, deltaTreshold = 1e-5, 
  nonnegative = false, maxIterations = 10000)

  N, M = size(X)
  lambda_eye = UniformScaling(lambda)
  xt = X'
  theta = xt * ((X * xt .+ lambda_eye) \ y)

  if nonnegative
    set_negative_to_zero!(theta)
  end

  iterations = 0

  while iterations < maxIterations
    iterations += 1

    # E-step:
    eta = theta

    # M-step:
    eta_squared = eta .^ 2
    xt_eta = broadcast(*, eta_squared, xt)
    theta = xt_eta * ((X * xt_eta .+ lambda_eye) \ y)

    # Ensure non-negativity if required
    if nonnegative
      set_negative_to_zero!(theta)
    end

    # We have converged if change is too small => break.
    if norm(theta - eta) < deltaTreshold
      break
    end
  end

  theta[abs(theta) .< epsilon] = 0.0

  theta
end

function l0EM_adaptus_solve(X, y, lambda = 2; 
  epsilon = 1e-4, deltaTreshold = 1e-5, 
  nonnegative = false, maxIterations = 10000)

  N, M = size(X)
  if N >= M
    lambda_eye = UniformScaling(lambda)
  else
    lambda_eye = lambda * eye(N)
  end
  xt = X'
  theta = xt * ((X * xt .+ lambda_eye) \ y)

  if nonnegative
    set_negative_to_zero!(theta)
  end

  iterations = 0

  while iterations < maxIterations
    iterations += 1

    # E-step:
    eta = theta

    # M-step:
    eta_squared = eta .^ 2
    xt_eta = broadcast(*, eta_squared, xt)
    theta = xt_eta * ((X * xt_eta .+ lambda_eye) \ y)

    # Ensure non-negativity if required
    if nonnegative
      set_negative_to_zero!(theta)
    end

    # We have converged if change is too small => break.
    if norm(theta - eta) < deltaTreshold
      break
    end
  end

  theta[abs(theta) .< epsilon] = 0.0

  theta
end

function set_negative_to_zero!(v)
  v[v .< 0.0] = 0.0
  v
end

function format_time(t)
  if t < 1e-11
    return @sprintf("%.2e nsec", t/1e-9)
  end
  sec_steps = [(1e-7, 1e-9, "nsec"), (1e-4, 1e-6, "usec"), 
    (1e-1, 1e-3, "msec"), (60, 1, "sec"), (3600, 60, "mins")]
  for i in 1:length(sec_steps)
    limit, divisor, name = sec_steps[i]
    if t < limit
      return @sprintf("%.2f %s", t/divisor, name)
    end
  end
  @sprintf("%.2f hours", t/3600)
end

using DataFrames
using HypothesisTests

coefficient_of_variation(v) = std(v) / mean(v)
mad(v) = median(abs(v .- median(v)))
rcv(v) = mad(v) / median(v)
iqr(v) = quantile(v, 0.75) - quantile(v, 0.25)
riqr(v) = iqr(v) / median(v)

SkipIndicator = -100
WP = symbol("p-value")
PM = symbol("+/- %")

summarize(df::AbstractDataFrame) = DataFrame(
  Avg = mean(df[:Time]), 
  PM = 100.0 * coefficient_of_variation(df[:Time]),
  #RCV = rcv(df[:Time]),
  #RIQR = riqr(df[:Time]),
  Max = maximum(df[:Time]),
  Min = minimum(df[:Time]),
  Reps = length(df[:Time]) 
)

#significance_indicator(p) = (p < 0.001) ? "***" : ((p < 0.01) ? "**" : ((p < 0.05) ? "*" : ((p == 1.0) ? "N/A" : "")))
significance_indicator(p) = (p < 0.001) ? "Yes (p < 0.001)" : ((p < 0.01) ? "Yes (p < 0.01)" : ((p < 0.05) ? "Yes (p < 0.05)" : ((p == 1.0) ? "N/A" : "No")))
speed_indication(t) = (t < 1) ? "." : ((t < 10) ? "," : ((t < 300) ? ":" : "|"))

calc_pvalue(fastest, times) = begin
  len = min(length(fastest), length(times))
  # We use a right sided p-value since we know times is at least as large as fastest, on average.
  pvalue(SignTest( (times .- fastest[1:len]), 0 ), tail = :right)
  # Wilcoxon/MannWhitneyUTest not really right to use since executions are paired...
  #pvalue(MannWhitneyUTest(times, fastest), tail = :both)
end

mean_speed_diff(fastest, times) = mean(times ./ fastest[1:length(times)])

summarize_measurements(times, indices, n, rep, descs, sortByAvg = false) = begin
  df = DataFrame({times, indices}, [:Time, :Function])

  # Skip rows that we have not yet filled in / measured:
  df = df[df[:Function] .> 0, :]

  # Skip rows that we have indicated as skipped
  df = df[!(df[:Time] .== SkipIndicator), :]

  # Summarize and calc relative speed
  dfs = by(df, :Function, summarize) 
  dfs[:Relative] = dfs[:Avg] / minimum(dfs[:Avg])

  # Find fastest function
  fastest = dfs[:Function][findfirst(dfs[:Avg], minimum(dfs[:Avg]))]

  # Compare all functions to the fastest with a wilcoxon test
  times_fastest = df[:Time][df[:Function] .== fastest]
  #dfs[WP] = map((j) -> pvalue(MannWhitneyUTest(times_fastest, df[:Time][df[:Function] .== j])), 1:n)
  dfs[WP] = map((j) -> calc_pvalue(times_fastest, df[:Time][df[:Function] .== j]), 1:n)
  
  # Calc relative for each paired rep and take mean of that:
  #dfs[:Rel2] = map((j) -> mean_speed_diff(times_fastest, df[:Time][df[:Function] .== j]), 1:n)

  # Make information more readable
  dfs[:Function] = descs
  if sortByAvg
    dfs = sort(dfs, cols = :Avg)
  end
  dfs[:Avg] = map(format_time, dfs[:Avg])
  #dfs[:Median] = map(format_time, dfs[:Median])
  dfs[symbol("Rel.")] = map((v) -> signif(v, 3), dfs[:Relative])
  map( (s) -> dfs[s] = map((v) -> signif(v, 2), dfs[s]), [WP] )
  dfs[:Max] = map(format_time, dfs[:Max])
  dfs[:Min] = map(format_time, dfs[:Min])
  ss = symbol("Signif slower?")
  dfs[ss] = map(significance_indicator, dfs[WP])
  dfs[PM] = map((v) -> signif(v, 3), dfs[:PM])

  dfs[:,[:Function, :Avg, symbol("Rel."), :Reps, PM, :Max, :Min, WP, ss]]
end

spaces(num) = join( map((i) -> " ", 1:num), "" )

# Compare the execution speed of a set of functions on the same randomly generated input.
# Run only as many replications of each function as is needed to show that it is slower than
# the fastest of the functions.
function speed_compare(funcs::Dict{ASCIIString, Function}, inputsGen::Function;
  MaxReplications = 15, 
  MinReplications = 3,
  alpha = 0.05)

  # Create description strings of equal length
  descs = collect(keys(funcs))
  max_desc_len = maximum(map(length, descs))
  formatted_descs = map((s) -> s * spaces(max_desc_len - length(s)), descs)

  fs = map((d) -> funcs[d], descs)
  n = length(fs)

  # Force JIT compilation
  inputs = inputsGen()
  map((i) -> fs[i](inputs...), 1:n)

  times = zeros(Float64, n * MaxReplications)
  indices = zeros(Int, n * MaxReplications)
  still_run = map((i) -> true, 1:n)

  local rep # since we need it in return statement below

  for rep in 1:MaxReplications
    inputs = inputsGen()
    for i in shuffle(collect(1:n))
      index = (i - 1) * MaxReplications + rep
      indices[index] = i
      if still_run[i]
        f = fs[i]
        times[index] = @elapsed f(inputs...)
        # println("  $(i). $(descs[i]): $(format_time(times[index]))")
      else
        times[index] = SkipIndicator
        # println("  $(i). $(descs[i]): skipped")
      end
      print(speed_indication(times[index]))
    end

    if rep >= MinReplications
      dfs = summarize_measurements(times, indices, n, rep, descs)
      # Only run those which do not have p-values below alpha compared
      # to the fastest running functions:
      still_run = still_run & (dfs[WP] .> alpha)
      if sum(still_run) < 2
        break
      end
    end
  end

  return summarize_measurements(times, indices, n, rep, formatted_descs, true)
end

alts = {
  "original" => l0EM_orig,
  "solve" => l0EM_solve,
  "unisc & solve" => l0EM_us_solve,
}
alts = convert(Dict{ASCIIString, Function}, alts)

input_gen_func(N, M, k, epsilon = 0.10) = () -> begin
  X = randn(N, M)
  #actual_theta = shuffle([linspace(1,k,k), zeros(M-k)][:])
  actual_theta = shuffle([10.0 * rand(k,1), zeros(M-k)][:])
  y = X * actual_theta + epsilon * randn(N, 1)
  (X, y)
end

alts2 = {
  "original" => l0EM_orig,
  "solve" => l0EM_solve,
  "unisc & solve" => l0EM_us_solve,
  # "ad. unisc & solve" => l0EM_adaptus_solve,
}
alts2 = convert(Dict{ASCIIString, Function}, alts2)

for N in [1e2, 1e3]
  for M in [1e2, 1e3, 1e4]
    for k in [4, 16]
      println("\nN = $(int(N)), M = $(int(M)), k = $(int(k))")
      showall( speed_compare(alts2, input_gen_func(int(N), int(M), k)) )
      println("")
    end
  end
end
