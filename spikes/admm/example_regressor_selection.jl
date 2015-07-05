# Cardinality constrained least-squares example (nonconvex)
include("regressor_selection_admm.jl")

function gen_sparse_reg_problem(m, n, p = 100/n)
    # generate sparse solution vector
    x = sprandn(n,1,p)

    # generate random data matrix
    A = randn(m,n)

    # normalize columns of A
    A = A ./ sqrt(sum(A.^2,1))

    # generate measurement b with noise
    b = A*x + sqrt(0.001)*randn(m,1)

    return A, x, b
end

# Ensure compiled
A, xtrue, b = gen_sparse_reg_problem(15, 50, 0.05);
xhat = regressor_sel(A, b, iceil(50*0.05), 1.0; MaxIterations = 10, Quiet = true);

mape(y, yhat) = mean(100.0 * abs((y .- yhat) ./ y))

function eval_sparse_reg(A, b, xtrue, xhat)
    bhat = A * xhat
    strue = Set(find(xtrue))
    shat = Set(find(xhat))
    missing = setdiff(strue, shat)
    common = setdiff(strue, missing)
    additional = setdiff(shat, strue)
    mapeval = mape(b, bhat)
    println("MAPE: $(round(mapeval,1))%\n",
        "  Support: common = $(length(common)), missing = $(length(missing)), additional = $(length(additional))")
    return [mapeval, length(common), length(missing), length(additional)]
end

srand(2)
m = 1500       # number of examples
n = 200       # number of features
p = 100/n      # sparsity density

A, xtrue, b = gen_sparse_reg_problem(m, n, p);
supportlen = length(find(xtrue));
#A, xtrue, b = map(fn -> readcsv(fn * ".csv"), ["A", "x", "b"]);
@time xhat = regressor_sel(A, b, p*n, 1.0; MaxIterations = 1000);
eval_sparse_reg(A, b, xtrue, xhat);
@time xhat2 = regressor_sel(A, b, supportlen, 1.0; MaxIterations = 1000);
eval_sparse_reg(A, b, xtrue, xhat2);

slens = collect(10 * (1:iceil(2*supportlen/10)));
numslens = length(slens)
res = zeros(4,numslens)
for i in 1:numslens
    xhat = regressor_sel(A, b, slens[i], 1.0; MaxIterations = 100);
    res[:,i] = map(float, eval_sparse_reg(A, b, xtrue, xhat))
end
