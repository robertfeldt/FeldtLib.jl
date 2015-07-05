# Cardinality constrained least-squares example (nonconvex)
include("regressor_selection.jl")

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
xhat = regressor_selection_admm(A, b, iceil(50*0.05), 1.0; MaxIterations = 10, Quiet = true);

include("eval_sparse_reg.jl")

srand(2)
m = 1500       # number of examples
n = 5000       # number of features
p = 100/n      # sparsity density

A, xtrue, b = gen_sparse_reg_problem(m, n, p);
supportlen = length(find(xtrue))
#A, xtrue, b = map(fn -> readcsv(fn * ".csv"), ["A", "x", "b"]);
@time xhat = regressor_selection_admm(A, b, p*n, 1.0; MaxIterations = 1000, Quiet = false);
eval_sparse_reg(A, b, xtrue, xhat);
@time xhat2 = regressor_selection_admm(A, b, supportlen, 1.0; MaxIterations = 1000, Quiet = false);
eval_sparse_reg(A, b, xtrue, xhat2);

slens = collect(10 * (1:iceil(2*supportlen/10)));
numslens = length(slens)
res = zeros(5,numslens)
for i in 1:numslens
    xhat = regressor_selection_admm(A, b, slens[i], 1.1; MaxIterations = 100);
    res[:,i] = map(float, eval_sparse_reg(A, b, xtrue, xhat))
end
