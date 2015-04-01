%function [X,obj,gap,iter,time] = smoothcovsel_solver(Sigma,rho,alpha,beta,maxiter,freqprint,eps)
% modified version to track progress
function [X,obj_vec,gap_vec,iter,time_vec] = smoothcovsel_solvermod(Sigma,rho,alpha,beta,maxiter,freqprint,eps)

obj_vec=[]; time_vec=[]; gap_vec=[];
start_t = cputime;
n = max(size(Sigma)); 
U_x = zeros(n, n);
U0 = U_x;

% initialize parameters 
[alpha,beta,L,N,sigma1] = initialize(U_x,Sigma,rho,alpha,beta,eps);
init_beta = beta;

% initialize sum_g %
sum_g = zeros(n,n);
total_k = 1; k = 0;

while (1==1)
    ttt=tic;
  % estimate upper bound on max eigenvalue 
  tmp = Sigma + rho*U_x; tmp = (tmp+tmp')*.5;
  [X, Lambda] = eig(tmp);
  vlambda = diag(Lambda);
  flag = 0;
  % safeguard against active iterates
  while (1==1)
    Lambda = zeros(n,1);
    for i=1:n
      if vlambda(i) <= 0
         Lambda(i) = beta;
      else
        Lambda(i) = min(max(1./vlambda(i), alpha), beta);
      end
    end
    lambdamin = min(Lambda); lambdamax = max(Lambda);
    if (lambdamax == beta) && (beta ~= init_beta) && (lambdamin > 0)
       beta = min(beta*1.05, init_beta); 
       flag = 1;
    else 
       break;
    end  
  end
  if flag == 1
    L = (rho*beta)^2; U0 = U_x;
    sum_g = zeros(n,n); k = 0;
  end
  % dynamic update on beta and parameters 
  if lambdamax/beta <= 0.95 
    beta = min(lambdamax*1.05, init_beta); 
    L = (rho*beta)^2; U0 = U_x; 
    sum_g = zeros(n,n); k = 0;
  end
  % compute gradient of phi(U) at U_x, and return X, fX, gap   
  % phi(U)=max{<-Sigma-rho*U, X> + logdet(X): alpha I <= X <= beta I} % 
  logdet = sum(log(Lambda));
  Lambda = diag(Lambda);
  X = X*Lambda*X';
  xsigma = sum(sum(X.*Sigma, 1));
  sumabsx = sum(sum(abs(X),1));
  fX = full(logdet - xsigma - rho*sumabsx);
  xu =  sum(sum(X.*U_x, 1));
  gap = full(rho*(sumabsx - xu));
  if mod(total_k,freqprint)==0
    disp(['Iter: ', num2str(total_k),'   Primal: ',num2str(fX,'%.10e'),'   Gap: ',num2str(gap,'%.10e')]);
  end
  if (gap <= eps) || (total_k >= maxiter)
    break;
  end
  % continue Nesterov's smooth optimization %
  g = -rho*X;
  U_y = max(min(U_x-g/L, 1), -1);% projection step%
  sum_g = sum_g + g*(k+1);
  U_x = max(min(U0-sum_g/(2*L), 1), -1);%prox step%
  U_x =(2*U_x+(k+1)*U_y)/(k+3);
  k = k+1; total_k = total_k+1;
  
  time_vec=[time_vec,toc(ttt)];
  obj_vec=[obj_vec,fX];
  gap_vec=[gap_vec,gap];
  
  
end
if (total_k >= maxiter)
  disp(['It exceeds the maximum iteration: ', num2str(maxiter)]);
end

% setup for output
obj = fX;
iter = total_k;
time = cputime-start_t;

% output iter and time %
outfile = fopen('out.txt', 'a');
fprintf(outfile, 'iter=%6d time=%7.1f obj=%e gap=%e\n',...
       [iter, time, obj, gap]);
fclose(outfile);
