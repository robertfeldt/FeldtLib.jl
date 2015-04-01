% ========================================================================
% Smooth optimization approach for one-norm penalized maximum
% likelihood problem:
% max  logdet(X) - <Sigma, X> - rho 1'|X|1 
% s.t. alpha I <= X <= beta I,
% where 
%    0 <= alpha < beta <= inf
% Note: beta can be inf; it also works for the case: alpha > 0, beta=inf. 
% ========================================================================
%
% [X,Obj,Gap,Iter,Time]=smoothcovsel_main(data_Sigma,rho,alpha,beta,maxiter,...
%                         freqprint,eps)
%
% Input:
% data_Sigma  - noisy covariance matrix
% rho         - penalty parameter
% alpha       - lower bound on eigenvalues 
% beta        - upper bound on eigenvalues
% maxiter     - maximum number of iterations
% freqprint   - printing frequency
% eps         - absolute accuracy of the solution
%
% Output:
% -------
%  X     - estimate of inverse of covariance matrix
%  Obj   - estimate of optimal function value
%  Gap   - optimality gap
%  Iter  - number of iterations
%  Time  - CPU time
%
function  [X, Obj, Gap, Iter, Time] = smoothcovsel_main(Sigma, rho,... 
                                        alpha, beta, maxiter, freqprint, eps)

% read the data matrix Sigma
%Sigma  - n x n matrix
%Sigma=load(data_Sigma);

[X, Obj, Gap, Iter, Time] = smoothcovsel_solver(Sigma, rho, alpha, beta,...
                               maxiter, freqprint, eps);

