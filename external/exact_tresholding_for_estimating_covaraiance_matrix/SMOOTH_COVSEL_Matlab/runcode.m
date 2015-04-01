% Test script for Covariance selection code.
% Initialize parameters ****************
% Data ....

ZERO=1.0e-4; % small perturbation
signoise=.15; % Signal to noise ratio
% Algorithm ....
rho=0.5;            % Controls sparsity
prec=1e-1;          % Numerical precision
numInst=10;          % Number of instances for tests

fid = fopen('result.txt','w');

for index=1:numInst

   n=100*index;         % Dimension
   nnz=(n*(n+1)*.5-n)*.01; % Number of nonzer coefficients in A
   % Form random test matrix 
   e=ones(n,n);
   rand('seed',20);
   % Generate A, the original inverse covariance, with random sparsity pattern...
   A=eye(n);
   for i=1:nnz
       A(ceil(n*rand),ceil(n*rand))=sign(rand()-.5);
   end
   A=A'*A;
   B=2*(rand(n,n)-.5*ones(n,n)); % Add noise...
   B=(B+B')*.5*signoise+inv(A+ZERO*eye(n));
   B=B-min([min(eig(B))-ZERO,0])*eye(n); % (Set min eig > 0 if noise level too high)
   B=(B+B')*.5;
  
   fprintf(fid,'---- test=%d n=%d rho=%3.1f ----\n',index,n,rho);
   % Test smooth optimization code with dynamic update on beta
   alpha = 0; beta = inf;
   [X,obj,gap,iter,time] = smoothcovsel_solver(B,rho,alpha,beta,inf,1,prec);
   fprintf(fid, 'iter=%d obj=%10.8f obj_gap=%10.8f time=%.2f\n', ...
            iter,obj,gap,time);
end

fclose(fid);

