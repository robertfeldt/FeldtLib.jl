# Julia version of Boyd et al's Matlab script available on:
#   https://web.stanford.edu/~boyd/papers/admm/nonconvex/regressor_sel.html

# regressor_sel  Solve lasso problem via ADMM
#
#  x, history = regressor_sel(A, b, K, rho, alpha)
#
# Attempts to solve the following problem via ADMM:
#
#   minimize || Ax - b ||_2^2
#   subject to card(x) <= K
#
# where card() is the number of nonzero entries.
#
# The solution is returned in the vector x.
#
# history is a structure that contains the objective value, the primal and
# dual residual norms, and the tolerances for the primal and dual residual
# norms at each iteration.
#
# rho is the augmented Lagrangian parameter.
#
# alpha is the over-relaxation parameter (typical values for alpha are
# between 1.0 and 1.8).
#
# More information can be found in the paper linked at:
# http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
#
function regressor_sel(A, b, K, rho;
    MaxIterations = 1000, Quiet = true, ABSTOL = 1e-4, RELTOL = 1e-2)

    # Data preprocessing and setup
    K = iceil(K)
    m, n = size(A)
    Atb = A' * b            # save a matrix-vector multiply

    local r_norm::Float64, s_norm::Float64, eps_pri::Float64, eps_dual::Float64
    x = zeros(n)
    z = zeros(n)
    u = zeros(n)
    L, U = factor(A, rho)   # cache the factorization

    for k in 1:MaxIterations
        # x-update
        q = Atb + rho*(z - u)    # temporary value
        if ( m >= n )    # if skinny
            x = U \ (L \ q)
        else            # if fat
            x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2
        end
        x = vec(x)

        # z-update with relaxation
        zold = z
        z = keep_largest!(x + u, K)

        # u-update
        u = u + (x - z)

        # Check convergence
        r_norm  = norm(x - z)
        s_norm  = norm(-rho*(z - zold))
        eps_pri = sqrt(n) * ABSTOL + RELTOL * max(norm(x), norm(-z))
        eps_dual= sqrt(n) * ABSTOL + RELTOL * norm(rho*u)

        if !Quiet
            println(@sprintf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n", k,
                 r_norm, eps_pri, s_norm, eps_dual, objective(A,b,x)))
        end

        if (r_norm < eps_pri) && (s_norm < eps_dual)
            break
        end
    end

    return z
end

function objective(A, b, x)
    sum( (A*x - b).^2 )
end

function keep_largest!(z, K)
    pos = sortperm(abs(z), rev=true)
    z[pos[K+1:end]] = 0
    z
end

import Base.factor
function factor(A, rho)
    m, n = size(A)
    if ( m >= n )  # if skinny
        L = chol( A'*A + rho*speye(n), :L )
    else            # if fat
        L = chol( speye(m) + 1/rho*(A*A'), :L )
    end

    return L, L'
end
