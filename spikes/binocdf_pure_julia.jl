# binocdf is the cdf for the Binomial distribution and binoccdf is its complement. 
# Distributions.jl implements the cdf of its Binomial by calling out to Rmath which 
# is a Julia package that wraps the corresponding functions used in R. Since we do not want to
# load the Distributions lib when starting Autotest we want a pure Julia implementation
# of binocdf. We base this code on the pure Ruby code for incomplete beta function available here:
#   http://rubydoc.info/gems/rubystats/0.2.3/Rubystats/SpecialMath#incomplete_beta-instance_method
#   

# Cumulative Distribution Function evaluated at x of the Binomial(n, p) distribution:
function binocdf(x, n, p)
  x = ifloor(x)
  if x == n
    return 1.0
  else
    incompletebeta(1-p, n - x, x + 1)
  end
end

# Complementary cdf of Binomial distr func:
binoccdf(x, n, p) = 1.0 - binocdf(x, n, p)

# incompletebeta(x, a, b) evaluates incomplete beta function, here a, b > 0 and 0 <= x <= 1.
function incompletebeta(x, a, b)
  if x == 0
    return 0.0
  elseif x == 1
    return 1.0
  else
    beta_gam = exp( a * log(x) + b * log(1.0 - x) - log_beta(a, b) )
    if x < ((a+1) / (a+b+2))
      return (beta_gam * beta_fraction(x, a, b) / a)
    else
      return (1.0 - (beta_gam * beta_fraction(1.0 - x, b, a) / b))
    end
  end
end

const LOG_GAMMA_X_MAX_VALUE = 2.55e292

function log_beta(a, b)
  if (a <= 0.0) || (b <= 0.0) || (a + b) > LOG_GAMMA_X_MAX_VALUE
    return 0.0
  else
    return (lgamma(a) + lgamma(b) - lgamma(a + b))
  end
end

const XMININ = 2.23e-303
const PRECISION = 8.88e-016

function beta_fraction(x, a, b, MaxIterations = 200) 
  c = 1.0
  sum_ab  = a + b
  a_plus  = a + 1.0
  a_minus = a - 1.0
  h = 1.0 - sum_ab * x / a_plus
  if abs(h) < XMININ
    h = XMININ
  end
  h     = 1.0 / h
  frac  = h
  m     = 1
  delta = 0.0

  while (m <= MaxIterations) && (abs(delta - 1.0) > PRECISION) 
    m2 = 2 * m
    # even index for d
    d = m * (b - m) * x / ( (a_minus + m2) * (a + m2))
    h = 1.0 + d * h
    if abs(h) < XMININ
      h = XMININ
    end
    h = 1.0 / h
    c = 1.0 + d / c
    if abs(c) < XMININ
      c = XMININ
    end
    frac *= (h * c)
    # odd index for d
    d = -(a + m) * (sum_ab + m) * x / ((a + m2) * (a_plus + m2))
    h = 1.0 + d * h
    if abs(h) < XMININ
      h = XMININ
    end
    h = 1.0 / h
    c = 1.0 + d / c
    if abs(c) < XMININ
      c = XMININ
    end
    delta = h * c
    frac *= delta
    m += 1
  end

  return frac
end

# Compare to Distributions.jl
using Distributions

# call once to compile things
binocdf(1, 10, 0.3)
binocdf(1.0, 10.0, 0.3)
binocdf(2.0, 10, 0.3)

# Sample lower values with higher prob through a log transform.
# Often more realistic since people tend to use lower ranges/values 
# frequently.
function logsample(min, max)
  shift = 1.0 - min
  logmin = log(min + shift)
  logmax = log(max + shift)
  exp(logmin + (logmax - logmin) * rand()) - shift
end

NumReps = 1e5
ourtimesum = 0.0
dtimesum = 0.0

for rep in 1:ifloor(NumReps)
  n = ifloor(logsample(2, 1e6))
  x = rand() * n # rand(0:n)
  p = rand()

  tic()
  ourres = binocdf(x, n, p)
  ourtimesum += toq()

  tic()
  dres = cdf(Binomial(n, p), x)
  dtimesum += toq()

  if !isapprox(ourres, dres)
    println("Failed for x = $x, n = $n, p = $p: our = $ourres, Distr.jl = $dres")
  end
end

println("Our average binocdf performance: ", ourtimesum / NumReps)
println("Distributions.jl average binocdf performance: ", dtimesum / NumReps)
println("Performance diff: ", round(ourtimesum / dtimesum, 2), "x")
