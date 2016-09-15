using Distributions

const AllContinuousDistributions = [
    Arcsine, Beta, BetaPrime, Cauchy, Chi, Chisq, Erlang, Exponential,
    FDist, Frechet, Gamma, GeneralizedExtremeValue, GeneralizedPareto,
    Gumbel, InverseGamma, InverseGaussian, Laplace, Levy, LogNormal,
    Logistic, Normal, NormalInverseGaussian, Pareto, Rayleigh, 
    SymTriangularDist, TDist, TriangularDist, Uniform, VonMises, Weibull
]

const ContinuousDistributionsThatCanFit = filter(AllContinuousDistributions) do D
    try
        fit_mle(D, rand(10))
        true
    catch _err
        false
    end
end

function fitalldistr{T<:Real}(x::AbstractArray{T}, sortby = :BIC)
    n = length(x)
    results = Any[]
    for D in ContinuousDistributionsThatCanFit
        try 
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
            push!(results, (d, criteria))
        catch _err
        end
    end
    sort(results, by = t -> t[2])
end

function fitall{T<:Real}(x::AbstractArray{T}, sortby = :BIC)
    fitalldistr(x, sortby)[1][1]
end
