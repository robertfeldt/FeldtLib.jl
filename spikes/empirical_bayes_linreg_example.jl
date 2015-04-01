# Example
#function ebreg_example(n=100, p=500, r=0.25, beta=0.6 * (1:5), M=2000)

n = 100
p = 500
r = 0.25
beta = 0.5 * (1:5)
numSamples = 2000
burnIn = 0.25
sampleCoefs = false

sigma2 = 1
alpha = 1
gamma = 0
logF = (x) -> dcomplex(x, n, p, 0.05, 1)

# Generate problem
s0 = length(beta)
beta0 = [beta, zeros(p-s0)]
S0 = float(beta0 .!= 0)
sig2_hat = nothing
X = randn(n, p)
#if(r != 0)
#  R = (1 - r) * diag(p) + array(r, c(p, p))
#  e = eigen(R)
#  sqR = e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
#  X = X * sqR
#end
y = X * beta0 + sqrt(sigma2) * randn(n)

# Call
@time models, sampledmodels = empirical_bayes_linreg(X, y; sigma2 = sigma2, alpha = alpha, gamma = gamma, 
  numSamples = int(1e4), logF = logF)

smodels = models[unique(sampledmodels)]

lengths = map((e) -> length(e[2]), models)
edges, counts = hist(float(lengths), 20)

logposteriors = map((e) -> e[6], models)
highestlogposterior = maximum(logposteriors)
bestmodelindex = find((m) -> m[6] == highestlogposterior, models)[1]
bestmodelentry = models[bestmodelindex]
bestmodel = bestmodelentry[3]

sses = map((e) -> e[5], models)
lowestsse = minimum(sses)
lowestssemodelindex = find((m) -> m[5] == lowestsse, models)[1]
lowestssemodelentry = models[lowestssemodelindex]

println("Number of sampled models: ", length(smodels))

println("Model with highest posterior: ", bestmodelentry)
println("\nModel with lowest sse: ", lowestssemodelentry)