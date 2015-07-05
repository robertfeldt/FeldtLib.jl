# Generalized ADMM function based on Boyd et al's Matlab scripts:
#   http://web.stanford.edu/~boyd/papers/admm/

abstract ADMMState
abstract NonconvexADMMState <: ADMMState

function admm_loop(s::ADMMState;
    MaxIterations = 1000,
    Quiet = true,
    ABSTOL::Float64 = 1e-4,
    RELTOL::Float64 = 1e-2)

    local r_norm::Float64, s_norm::Float64, eps_pri::Float64, eps_dual::Float64
    local best_objval::Float64 = Inf # Make sure its updated in 1st iteration
    local best_solution

    for k in 1:MaxIterations
        zold = s.z
        update_xzu!(s, k)

        # Check convergence
        r_norm  = norm(s.x - s.z)
        s_norm  = norm(-s.rho*(s.z - zold))
        eps_pri = sqrt(s.n) * ABSTOL + RELTOL * max(norm(s.x), norm(-s.z))
        eps_dual= sqrt(s.n) * ABSTOL + RELTOL * norm(s.rho*s.u)
        objval = objective(s)

        if !Quiet
            println(@sprintf("%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f", k,
                 r_norm, eps_pri, s_norm, eps_dual, objval))
        end

        if objval < best_objval
            best_objval = objval
            best_solution = copy(solution(s))
        end

        if (r_norm < eps_pri) && (s_norm < eps_dual)
            break
        end
    end

    # Return the solution with best objective value if this is non-convex ADMM.
    if typeof(s) <: NonconvexADMMState
        return best_solution
    else
        return solution(s)
    end
end
