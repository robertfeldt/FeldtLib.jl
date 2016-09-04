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

function optim_based_slope_regression(y, X; q = 0.01, sigma = 0.1, absmax = nothing)
    n, p = size(X)
    @assert n == length(y)
    lambdas = lambdas_bh(p, q)

    if typeof(absmax) <: Number
        # Can't use Fminbox in Optim when also have Stan.jl loaded since they both export Newton...
        #lower = (-absmax) * ones(p)
        #upper = abs(lower)
        #res = optimize(f, randn(p)*sigma, lower, upper, Fminbox(), optimizer = LBFGS(), 
        #        OptimizationOptions(iterations = 5000))

        # So doing an explicit extra pressure to not go outside the absmax box
        function f(betas)
            objval = slope_objective_func(betas, y, X, lambdas, sigma)
            abetas = abs(betas)
            outside_box = abetas .> absmax
            if length(outside_box) > 0
                objval += sumabs2(betas[outside_box])
            end
            objval
        end
        res = optimize(f, randn(p)*sigma, SimulatedAnnealing(), 
                OptimizationOptions(iterations = 20000))
    else
        f(betas) = slope_objective_func(betas, y, X, lambdas, sigma)
        res = optimize(f, randn(p)*sigma, LBFGS(), 
                OptimizationOptions(iterations = 5000))
    end

    return Optim.minimizer(res) 
end