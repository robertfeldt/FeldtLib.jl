# The normal fit api in Distributions.jl is:
#   fit(Distr, values)
# so we want to do a more general
#   fitall(values)
# that will try all distributions that can be fitted and return back the "best" one
# for relevant criteria of "best", typically BIC or AIC.

using Distributions

AllContinuousDistributions = [
    Arcsine, Beta, BetaPrime, Cauchy, Chi, Chisq, Erlang, Exponential,
    FDist, Frechet, Gamma, GeneralizedExtremeValue, GeneralizedPareto,
    Gumbel, InverseGamma, InverseGaussian, Laplace, Levy, LogNormal,
    Logistic, Normal, NormalInverseGaussian, Pareto, Rayleigh, 
    SymTriangularDist, TDist, TriangularDist, Uniform, VonMises, Weibull
]

samples1 = rand(Normal(3.4, 1.2), 100)

ContinuousDistributionsThatCanFit = filter(AllContinuousDistributions) do D
    try
        fit_mle(D, samples1)
        true
    catch _err
        false
    end
end

# Start for continuous distributions
function fitalldistr{T<:Real}(x::AbstractArray{T}, sortby = :BIC)
    n = length(x)
    allfits = map(ContinuousDistributionsThatCanFit) do D
        d = fit_mle(D, x)
        loglik = loglikelihood(d, x)
        k = length(params(d))
        if sortby == :BIC
            criteria = -2 * loglik + k * log(n)
        elseif sortby == :AIC
            criteria = -2 * loglik + 2 * k
        elseif sortby == :AICc
            criteria = -2 * loglik + 2 * k + ((2*k*(k+1))/(n-k-1))
        elseif sortby == :NLogL
            criteria = -loglik
        end
        (d, criteria)
    end
    sort(allfits, by = t -> t[2])
end

function fitall{T<:Real}(x::AbstractArray{T}, sortby = :BIC)
    fitalldistr(x, sortby)[1][1]
end

samples1 = rand(Normal(3.4, 1.2), 100)
d1 = fitall(samples1)

samples2 = rand(Exponential(0.3), 100)
d2 = fitall(samples2)

samples3 = rand(LogNormal(-5.0, 2.0), 100)
d3s = fitalldistr(samples3)
d3 = fitall(samples3)

samples4 = rand(Gamma(6.20, 1.3), 200)
d4s = fitalldistr(samples4)
d4 = fitall(samples4)

# Should we maybe use a non-parametric test to select between them?