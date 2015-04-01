function [Time, X, S,C] = compute_GL_exact(Sigma_cov,rho,MAX_COMP_SIZE,maxiter,freqprint,epstol)
%[Time, X, S,C] = compute_GL_exact(Sigma_cov,rho,MAX_COMP_SIZE,maxiter,freqprint,epstol)
%
%This function performs the exact thresholding of covariance-matrices in
%the context of the Graphical Lasso problem:
% \mini_{X psd}  -\log\det(X) + \trace(Sigma_cov X) + \rho *
% sum(sum(abs(X)))  ----- (GL)
%
% For a large value of rho, the solution X* to (GL)
% often breaks up into connected components. These components can be exactly 
% identified by simple operations on the sample covariance matrix i.e. Sigma_cov
%
% INPUTS: 
% Sigma_cov : Sample covariance matrix (p \times p matrix), PSD 
% rho       : penalty parameter
%MAX_COMP_SIZE: Maximal allowable component-size (dependent upon rho)
%               defaults to 1K    
% maxiter,freqprint,epstol : these are all arguments for
% `smoothcovsel_main' ( solver for (GL), By Z. Lu  
% downloaded from http://people.math.sfu.ca/~zhaosong/ )
%   maxiter: max number of iters to perform (default: 500)
%   freqprint: (integer) to print progress in after as many iters (default: maxiter + 1)
%   epstol: convergence criterion. (default: 1e-4)
%[for details please see the documentation of the solver in folder
% /SMOOTH_COVSEL_Matlab]
%
% OUTPUTS:
% X : precision matrix i.e. the inverse covariance matrix
% Time: Times in secs taken for each of the connected components
% S : Number of connected components in solution
% C : labels of component memberships

% code written by Rahul Mazumder
%For any questions/ suggestions/ comments/ bugs please report to rahulm@stanford.edu


if nargin < 2
 disp('number of input arguments must be at least 2');
end

if (isempty(Sigma_cov))
    disp('Sigma_cov must be provided');
end

if (isempty(rho))
disp('rho must be provided');
end
    
if (nargin == 2) 
MAX_COMP_SIZE=[];maxiter=[];freqprint=[];epstol=[];
elseif (nargin ==3) 
 maxiter=[];freqprint=[];epstol=[];
elseif (nargin == 4) 
freqprint=[];epstol=[];
else
epstol=[];
end

 
if isempty(maxiter)
maxiter=500; 
end


if (isempty(freqprint))
freqprint=maxiter+1;
end

if (isempty(epstol))
    epstol=10^-4;
end


if (size(Sigma_cov,1) ~= size(Sigma_cov,2))
disp('Sigma_cov needs to be square')
end



[Time, X,S,C] = compute_GL_exact_smacs(Sigma_cov,rho,MAX_COMP_SIZE,maxiter,freqprint,epstol);


end

    
