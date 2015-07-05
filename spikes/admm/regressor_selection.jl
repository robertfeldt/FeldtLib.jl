include("admm_common.jl")

type RegressorSelectionState <: NonconvexADMMState
    rho::Float64
    K::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    m::Int
    n::Int
    At::Matrix{Float64}
    Atb::Vector{Float64}
    x::Vector{Float64}
    z::Vector{Float64}
    u::Vector{Float64}
    L::Matrix{Float64}
    U::Matrix{Float64}
    function RegressorSelectionState(A::Matrix{Float64}, b::Vector{Float64},
        K::Int, rho::Float64)
        m, n = size(A)
        At = A'
        Atb = At * b            # save a matrix-vector multiply
        L, U = factor(A, rho)   # cache the factorization
        new(rho, iround(K), A, b, m, n, At, Atb,
            zeros(n), zeros(n), zeros(n), L, U)
    end
end

function regressor_selection_admm(A, b, K, rho::Float64 = 1.0; options...)
    as = RegressorSelectionState(A, vec(b), iround(K), rho)
    admm_loop(as; options...)
end

function update_xzu!(s::RegressorSelectionState, iter::Int)
    # x-update
    q = s.Atb + s.rho*(s.z - s.u)    # temporary value
    if ( s.m >= s.n )    # if skinny
        x = s.U \ (s.L \ q)
    else            # if fat
        x = q/s.rho - (s.At*(s.U \ ( s.L \ (s.A*q) )))/s.rho^2
    end
    s.x = vec(x)

    # z-update with relaxation
    s.z = keep_largest!(s.x + s.u, s.K)

    # u-update
    s.u = s.u + (s.x - s.z)
end

solution(s::RegressorSelectionState) = s.z

function objective(s::RegressorSelectionState)
    sum( (s.A*s.x - s.b).^2 )
end

function keep_largest!(z::Vector{Float64}, K::Int)
    pos = sortperm(abs(z), rev=true)
    z[pos[(K+1):end]] = 0
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
