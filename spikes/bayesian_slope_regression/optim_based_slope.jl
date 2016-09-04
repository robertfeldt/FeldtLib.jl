using Optim

# SLOPE objective to be minimized
function slope_objective_func(betas, y, X, lambdas, sigma)
    0.5 * sumabs2(y .- X*betas) + sigma * sum(sort(betas, rev=true) .* lambdas)
end

# Set lambdas according to Benjamini-Hochberg
using Distributions
const NormalDistr = Normal(0, 1)
function lambdas_bh(p, q)
    lambdas = zeros(Float64, p)
    for i in 1:p
        lambdas[p+1-i] = quantile(NormalDistr, 1-0.5*i*q/p)
    end
    lambdas
end

function optim_based_slope_regression(y, X, q = 0.01, sigma = 0.1)
    n, p = size(X)
    @assert n == length(y)
    lambdas = lambdas_bh(p, q)
    f(betas) = slope_objective_func(betas, y, X, lambdas, sigma)
    #res = optimize(f, zeros(p), NelderMead())
    res = optimize(f, randn(p)*sigma, LBFGS(), 
            OptimizationOptions(iterations = 5000))
    return Optim.minimizer(res) 
end