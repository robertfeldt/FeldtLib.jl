using Mamba, Stan

currdir = (@__FILE__() == nothing) ? "." : dirname(@__FILE__())
const BayesianSlopeModel = open(fh -> readall(fh), joinpath(currdir, "BayesianSlope.stan"))

# Create fake data to test
NumCoeffs = 100
nonzero_coeffs = [2, 3, -4, 7]
nx = length(nonzero_coeffs)
NumSpurious = NumCoeffs - nx
coeffs = zeros(NumCoeffs)
idxs_nonzero = shuffle(collect(1:NumCoeffs))[1:nx]
coeffs[idxs_nonzero] = nonzero_coeffs
p = length(coeffs)
s = 0.20
n = 200
X = randn(n, p)
y = X * coeffs .+ s * randn(n)

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

# Prep and then call stan
data = Dict("p" => p, "n" => n, "y" => y, "X" => X, 
    "lambda" => lambdas_bh(p, 0.01), "a" => 1.0, "g" => 1.0);
bayes_slope_data = [data];

num_betas = length(coeffs);
beta_names = map(i -> "beta.$i", 1:num_betas);
convergence_stats = ["lp__", "accept_stat__"];

stanmodel = Stanmodel(name="bayesian_slope", model=BayesianSlopeModel, 
    monitors = vcat(convergence_stats, beta_names), update=2000)
@time sim1 = stan(stanmodel, bayes_slope_data, currdir, CmdStanDir=CMDSTAN_HOME);

convstats_chain1 = sim1[1:size(sim1, 1), convergence_stats, 1];

betas_chain1 = sim1[1:size(sim1, 1), beta_names, 1].value;
est_betas_chain1 = median(betas_chain1, 1)[:];
#std_betas_chain1 = std(betas_chain1, 1)[:];
#includes_zero(val, stdest) = (val - stdest) < 0.0 < (val + stdest)
#nonzero_idxs = filter(i -> !includes_zero(est_betas_chain1[i], std_betas_chain1[i]), 1:length(est_betas_chain1));
includes_zeroes_hdr(vals, level=0.05) = (qs = quantile(vals, [level, 1.0-level]); qs[1] < 0.0 < qs[2])
nonzero_idxs = filter(i -> !includes_zeroes_hdr(betas_chain1[:,i,:][:]), 1:length(est_betas_chain1));
println("Non-zero coeffs:")
map(i -> println("  $i: ", join(map(v -> round(v,2), quantile(betas_chain1[:,i,:][:], [0.05, 0.50, 0.95])), " - ")), sort(nonzero_idxs));
println("Actual coeffs:")
map(i -> println("  $i: ", coeffs[i]), sort(idxs_nonzero));

# Now compare to Optim-based SLOPE regression
include("optim_based_slope.jl")
betas_optim = optim_based_slope_regression(y, X)
largest_abs = sortperm(abs(betas_optim), rev=true)
println("Largest 10: ", largest_abs[1:10])
println("Coeffs:     ", map(i -> round(betas_optim[i], 3), largest_abs[1:10]))

